import numpy as np
import pandas as pd
import random
import math
import os

from blood import *

class Request:

    def __init__(self, day_needed, day_available, num_units, patgroup, blood):
        self.day_needed = day_needed
        self.day_available = day_available
        self.num_units = num_units
        self.patgroup = patgroup
        self.blood = blood

class Demand():
    """
    Mean and Coefficient of Variation (CV) for each day of the week.
    """
    def __init__(self, SETTINGS, weekday_index):

        daily_demand_distribution = np.array(
            [[2.160824e+03, 7.144959e-02],
            [1.940279e+03,  8.390841e-02],
            [1.980793e+03,  7.903061e-02],
            [1.954597e+03,  7.880648e-02],
            [2.029441e+03,  7.298936e-02],
            [319.5225256,   0.2143665],
            [308.0818815,   0.2747322]])

        avg_weekly_demand = np.sum(daily_demand_distribution, axis=0)[0]

        assert (weekday_index >= 0 and weekday_index <= 6), "Day index must be between 0 and 6"

        # number of identical subdistributions which together form exactly the daily demand distribution
        sizefactor = avg_weekly_demand / SETTINGS.avg_daily_demand

        national_mean = daily_demand_distribution[weekday_index, 0]
        national_cv = daily_demand_distribution[weekday_index, 1]
        national_stdev = national_cv * national_mean

        # compute mean and stdev for subdistribution (= scaled down distribution)
        mean = national_mean / sizefactor
        stdev = national_stdev / math.sqrt(sizefactor)

        # procedure below is based on Source: https://www.win.tue.nl/~iadan/alqt/fit.pdf
        # compute parameters for fitting
        self.a = (stdev / mean)**2 - (1 / mean)

        # mixed geometric distribution
        if self.a >= 1:
            self.p1 = (mean * (1 + self.a + math.sqrt(self.a**2 - 1))) / (2 + mean * (1 + self.a + math.sqrt(self.a**2 - 1)))
            self.p2 = (mean * (1 + self.a - math.sqrt(self.a**2 - 1))) / (2 + mean * (1 + self.a - math.sqrt(self.a**2 - 1)))
            self.q1 = 1 / (1 + self.a + math.sqrt(self.a**2 - 1))
        else:
            self.k = math.floor(1 / self.a)
            self.q = ((self.k + 1) * self.a - math.sqrt((self.k + 1) * (1 - self.a * self.k))) / (1 + self.a)
            self.p = mean / ((self.k + 1) - self.q + mean)


        self.patgroups = ["Other", "Wu45", "MDS", "Thal", "AIHA", "ALA", "SCD"]

        if SETTINGS.demand_scenario == "Other":
            self.patgroup_distr = [1, 0, 0, 0, 0, 0, 0]
        elif SETTINGS.demand_scenario == "regional":
            self.patgroup_distr = [0.886867, 0.049352, 0, 0.008607969, 0.014933, 0.031632, 0.008607969]
        elif SETTINGS.demand_scenario == "university":
            self.patgroup_distr = [0.64605, 0.10250, 0.04542, 0.05665, 0.02731, 0.06543, 0.05665]


    def sample_number_of_units(self):
        """
        Sample a number of units demanded based on parameters as computed in the constructor of this class
        Source: https://www.win.tue.nl/~iadan/alqt/fit.pdf
        For sampling the geometric distribution we use (1-p) as parameter and subtract one to account for the difference in 
        definitions for the geometric distribution. 
        See https://en.wikipedia.org/wiki/Geometric_distribution#:~:text=The%20geometric%20distribution%20gives%20the,%2C%203%2C%20
        Paper uses the definition on the right, whereas the implementation uses the definition on the left.
        """
        if self.a >= 1: #sample geometric
            if random.random() < self.q1:
                return self.sample_geometric(1 - self.p1) - 1
            else:
                return self.sample_geometric(1 - self.p2) - 1
        else: #sample negative binomial (repeated geometric sampling)
            num_units_demanded = 0
            if random.random() < self.q:
                i = 0
                while i < self.k:
                    num_units_demanded += self.sample_geometric(1 - self.p) - 1
                    i += 1
            else:
                i = 0
                while i < self.k+1:
                    num_units_demanded += self.sample_geometric(1 - self.p) - 1
                    i += 1
            
            return num_units_demanded

    # Generate a list of random requests according to the distribution which was calculated in the constructor of this instance
    def sample_request_for_day(self, SETTINGS, day_index, df):
        num_units_demanded = self.sample_number_of_units()

        num_units = 0
        while num_units < num_units_demanded:
            r = self.get_random_request(SETTINGS, day_index)
            df.loc[len(df)] = [r.day_needed, r.day_available, r.num_units, r.patgroup, r.blood.ethnicity] + r.blood.vector
            num_units += r.num_units

        return df

    def get_random_request(self, SETTINGS, day_needed):

        patgroup = random.choices(self.patgroups, weights = self.patgroup_distr, k=1)[0]
        lead_time = random.choices(range(14), weights = SETTINGS.request_lead_time_probabilities[patgroup], k=1)[0]
        num_units = random.choices(range(1,5), weights = SETTINGS.request_num_units_probabilities[patgroup], k=1)[0]

        if patgroup == "SCD":
            ethnicity = "African"
        else:
            ethnicity = "Caucasian"

        return Request(day_needed, max(0, day_needed - lead_time), num_units, patgroup, Blood(SETTINGS, ethnicity))


    # Sample a geometric distribution with parameter p (mean = 1/p)
    def sample_geometric(self, p):
        return math.ceil(math.log(1 - random.random()) / math.log(1 - p))



# Generate multiple demand scenarios with the same parameters.
def generate_demand(SETTINGS):

    avg_daily_demand = SETTINGS.avg_daily_demand
    duration = SETTINGS.test_days + SETTINGS.init_days
    name = SETTINGS.name

    path = SETTINGS.home_dir + "demand"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"demand/{avg_daily_demand}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    i = 0
    while os.path.exists(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{name}_{i}.csv"):
        i += 1

    for _ in range(SETTINGS.episodes):

        print(f"Generating demand '{name}_{i}'.")

        # initialize distributions for each day of the week
        daily_distributions = []
        for weekday_index in range(7):
            daily_distributions.append(Demand(SETTINGS, weekday_index))

        df = pd.DataFrame(columns = ["Day Needed", "Day Available", "Num Units", "Patient Type", "Ethnicity"] + SETTINGS.major + SETTINGS.minor)

        for day_index in range(duration):
            df = daily_distributions[day_index%7].sample_request_for_day(SETTINGS, day_index, df)

        df.to_csv(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{name}_{i}.csv", index=False)

        i += 1





#     class DemandScenario:

#         #/ <summary>
#         #/ Create Demand scenario by reading it from csv file
#         #/ </summary>
#         #/ <param name="duration"></param>
#         #/ <param name="avg_daily_demand"></param>
#         #/ <param name="name"></param>
#         #/ <param name="normal_patients_one_day_ahead"></param>
#         #/ <param name="remove_twodotfivesup"></param>
#         # public DemandScenario(int duration, int avg_daily_demand, string name, bool normal_patients_one_day_ahead = false, bool remove_twodotfivesup = true)
#         def __init__(self, duration, avg_daily_demand, name):
#             #instance fields found by C# to Python Converter:
#             self.Duration = 0
#             self.Name = None
#             self.Requests = None
#             self.indicesPerDay = []
#             self.AverageDailyDemand = 0
#             self.PatientEthnicityDist = []

#             self.Duration = duration
#             self.Name = name
#             self.AverageDailyDemand = avg_daily_demand
#             try:
#                 # this.Requests = File.ReadAllLines("../../../Multi-day-demand-scenarios/" + avg_daily_demand + "/" + Duration + "/" + name + ".csv").Skip(1).Select(v => Request.FromCsv(v)).ToList()
#                 self.Requests = File.ReadAllLines("../../../demand/" + str(avg_daily_demand) + "/" + str(self.Duration) + "/" + name + ".csv").Skip(1).Select(lambda v : Request.FromCsv(v)).ToList()
#             except Exception as e:
#                 print("An exception occured when loading demand scenario " + name)
#                 raise e

#             #set ethnicity dist
#             self.PatientEthnicityDist = [0 for _ in range(3)]
#             totalNumberOfPatients = len(self.Requests)
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#             self.PatientEthnicityDist[0] = self.Requests.Where(lambda r : r.Blood.Ethnicity == Ethnicity.Caucasian).Count() / float(totalNumberOfPatients)
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#             self.PatientEthnicityDist[1] = self.Requests.Where(lambda r : r.Blood.Ethnicity == Ethnicity.Negroid).Count() / float(totalNumberOfPatients)
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#             self.PatientEthnicityDist[2] = self.Requests.Where(lambda r : r.Blood.Ethnicity == Ethnicity.Asian).Count() / float(totalNumberOfPatients)

#             for r in self.Requests:
#                 r.Blood.PrecomputeUsabilities(self)

#             #precompute indicesPerDay
#             self.indicesPerDay = [None for _ in range(duration)]
#             for i in range(0, duration):
#                 self.indicesPerDay[i] = []

#             r = 0
#             while r < len(self.Requests):
#                 self.Requests[r].Index = r
#                 t = Requests[r].DayAvailable
#                 while t <= self.Requests[r].DayNeeded:
#                     self.indicesPerDay[t].append(r) #add index of request r to each day at which the request is known
#                     t += 1
#                 # if ((int)Requests[r].PatientType<=1 && Requests[r].DayNeeded - 1 >= 0 && normal_patients_one_day_ahead)
#                 #     this.indicesPerDay[Requests[r].DayNeeded - 1].Add(r)
#                 r += 1

#             # for (int i = 0; i < duration; i++)
#             # {
#             #     if (remove_twodotfivesup && (float)indicesPerDay[i].Select(index => Requests[index]).Where(r => r.DayNeeded == i).Sum(r => r.NumUnits) > avg_daily_demand * 2.5f)
#             #     {
#             #         List<Request> today = indicesPerDay[i].Select(a => Requests[a]).Where(r=> r.DayNeeded == i).ToList()
#             #         int sum = today.Sum(r => r.NumUnits)
#             #         List<Request> newlist_1 = indicesPerDay[i].Where(r_index => Requests[r_index].DayNeeded != i).Select(x => Requests[x]).ToList()
#             #         List<int> newlist_2 = newlist_1.Select(x => x.Index).ToList()
#             #         List<int> curday = indicesPerDay[i].Where(r_index => Requests[r_index].DayNeeded == i).ToList()
#             #         float total = 0
#             #         int index = 0
#             #         while (total < avg_daily_demand * 2.5f)
#             #         {
#             #             newlist_2.Add(curday[index])
#             #             total += Requests[curday[index]].NumUnits
#             #             index++
#             #         }
#             #         indicesPerDay[i] = newlist_2
#             #     }

#             # }


#         #/ <summary>
#         #/ Get all requests available at a given day
#         #/ </summary>
#         #/ <param name="day_index"></param>
#         #/ <returns></returns>
#         def GetDemandForDay(self, day_index):
#             indices = self.indicesPerDay[day_index]
#             return indices.Select(lambda r : self.Requests[r]).ToList()

#         #/ <summary>
#         #/ Count the total number of requests that occur during the given duration
#         #/ </summary>
#         #/ <param name="duration"></param>
#         #/ <returns></returns>
#         def GetTotalNumberOfRequests(self, duration, init_length):
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#             return self.Requests.Where(lambda r : r.DayNeeded < duration and r.DayNeeded>= init_length).Count()

#         #/ <summary>
#         #/ Limit the demand of each day to maximum_relative_demand * AverageDailyDemand units
#         #/ </summary>
#         #/ <param name="maximum_relative_demand">size of demand relative to AverageDailyDemand. value 1 = AverageDailyDemand</param>
#         def RemoveToLargeDemandDays(self, maximum_relative_demand):
#             #loop over all days
#             i = 0
#             while i < self.Duration:
#                 #if total number of units requested on current day is larger than the allowed maximum
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#                 if float(self.indicesPerDay[i].Select(lambda index : self.Requests[index]).Where(lambda r : r.DayNeeded == i).Sum(lambda r : r.NumUnits)) > self.AverageDailyDemand * maximum_relative_demand:
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#                     today = self.indicesPerDay[i].Select(lambda a : self.Requests[a]).Where(lambda r : r.DayNeeded == i).ToList()
#                     sum = today.Sum(lambda r : r.NumUnits)
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#                     newlist_1 = self.indicesPerDay[i].Where(lambda r_index : self.Requests[r_index].DayNeeded != i).Select(lambda x : self.Requests[x]).ToList()
#                     newlist_2 = newlist_1.Select(lambda x : x.Index).ToList()
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to LINQ queries:
#                     curday = self.indicesPerDay[i].Where(lambda r_index : self.Requests[r_index].DayNeeded == i).ToList()

#                     #construct new list of requests and only allow requests if the current total is lower than the allowed maximum.
#                     total = 0
#                     index = 0
#                     while total < self.AverageDailyDemand * maximum_relative_demand:
#                         newlist_2.append(curday[index])
#                         total += self.Requests[curday[index]].NumUnits
#                         index += 1
#                     self.indicesPerDay[i] = newlist_2
#                 i += 1

#         #/ <summary>
#         #/ Create a Demand Scenario and write it to file
#         #/ </summary>
#         #/ <param name="average_daily_demand">The average daily demand (= week average divided by 7) </param>
#         #/ <param name="name">Name of the scenario for storing</param>
#         #/ <param name="duration">Duration in days, by default 396=365+31</param>
#         @staticmethod
#         def CreateDemandScenarioAndWriteToCSV(average_daily_demand, name, patientEthnicityDist, duration = 396):
#             #initialize distributions for each day of the week
#             daily_distributions = [None for _ in range(7)]
#             for day_index in range(0, 7):
#                 daily_distributions[day_index] = DemandDistribution(day_index, average_daily_demand)

#             #sample requests for every day
#             requests = []
#             for i in range(0, duration):
#                 requests.extend(daily_distributions[math.fmod(i, 7)].SampleRequestsForDay(i, patientEthnicityDist))

#             print("Creating Demand Scenario '" + name + "' with average daily demand " + str(average_daily_demand) + " and duration " + str(duration))

#             #write all requests to csv file
#             with System.IO.StreamWriter("../../../demand/" + str(average_daily_demand) + "/" + str(duration) + "/" + name + ".csv") as file:
#                 header = ["Day Needed", "Day Available", "Num Units", "Patient Type", "Ethnicity"]
#                 header.extend(Blood.Antigens)
#                 file.WriteLine(string.Join(",", header.ToArray()))
#                 file.Flush()

#                 for r in requests:
#                     file.WriteLine(Request.ToCSVString(r))


#         