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
    def __init__(self, SETTINGS, weekday_index, avg_daily_demand):

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
        sizefactor = avg_weekly_demand / avg_daily_demand

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
    def sample_request_for_day(self, SETTINGS, PARAMS, day_index, df, htype):
        num_units_demanded = self.sample_number_of_units()

        num_units = 0
        while num_units < num_units_demanded:
            r = self.get_random_request(SETTINGS, PARAMS, day_index, htype)
            df.loc[len(df)] = [r.day_needed, r.day_available, r.num_units, r.patgroup, r.blood.ethnicity] + r.blood.vector
            num_units += r.num_units

        return df

    def get_random_request(self, SETTINGS, PARAMS, day_needed, htype):

        patgroup = random.choices(PARAMS.patgroups, weights = [PARAMS.patgroup_distr[htype][pg] for pg in PARAMS.patgroups], k=1)[0]
        lead_time = random.choices(range(14), weights = PARAMS.request_lead_time_probabilities[patgroup], k=1)[0]
        num_units = random.choices(range(1,5), weights = PARAMS.request_num_units_probabilities[patgroup], k=1)[0]

        if patgroup == "SCD":
            ethnicity = "African"
        else:
            ethnicity = "Caucasian"

        return Request(day_needed, max(0, day_needed - lead_time), num_units, patgroup, Blood(PARAMS, ethnicity))


    # Sample a geometric distribution with parameter p (mean = 1/p)
    def sample_geometric(self, p):
        return math.ceil(math.log(1 - random.random()) / math.log(1 - p))



# Generate multiple demand scenarios with the same parameters.
def generate_demand(SETTINGS, PARAMS, htype, avg_daily_demand):

    duration = SETTINGS.test_days + SETTINGS.init_days

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
    while os.path.exists(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{htype}_{i}.csv"):
        i += 1

    for _ in range(SETTINGS.episodes[0],SETTINGS.episodes[1]):

        print(f"Generating demand '{htype}_{i}'.")

        # initialize distributions for each day of the week
        daily_distributions = []
        for weekday_index in range(7):
            daily_distributions.append(Demand(SETTINGS, weekday_index, avg_daily_demand))

        df = pd.DataFrame(columns = ["Day Needed", "Day Available", "Num Units", "Patient Type", "Ethnicity"] + PARAMS.major + PARAMS.minor)

        for day_index in range(duration):
            df = daily_distributions[day_index%7].sample_request_for_day(SETTINGS, PARAMS, day_index, df, htype)

        df.to_csv(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{htype}_{i}.csv", index=False)

        i += 1



    # class PatientType(Enum):
    #     OTHER = 0
    #     WU45 = 1
    #     MDS = 2
    #     THAL = 3
    #     AIHA = 4
    #     ALA = 5      
    #     SCD = 6

    # class PatientTypeDistributionOrigin(Enum):
    #     ALLNORMAL = 0
    #     OLVG = 1
    #     AMC = 2

    


        #variables used for multi-hospital setup

#         # -1 is must
#         PatientTypeCostMatrix = [[-1, -1, -1, .02651, .05429, .18433, .06439, .29543, 0, 0, 0, .00337, .00128, .02273, .00673, .04293, .00168], 
#         [-1, -1, -1, .02651, -1, -1, .06439, -1, 0, 0, 0, .00337, .00128, .02273, .00673, .04293, .00168], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .10100, .03844, .68177, .20200, 1.28778, .05050], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .10100, .03844, .68177, .20200, 1.28778, .05050], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .20200, .07689, 1.3653, .40401, 2.57556, .10100], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .20200, .07689, 1.3653, .40401, 2.57556, .10100], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .33667, .12815, -1, .67335, -1, -1]]

        

#         #each of these arrays specifies the probability that a patient is of a specific category (in order of the PatientType enum above)
#         PatientTypeProbabilitiesOLVG = [.886867, .049352, 0, .008607969, .014933, .031632, .008607969]
#         PatientTypeProbabilitiesAMC = [.64605, .10250, .04542, .05665, .02731, .06543, .05665]
#         PatientTypeProbabilitiesAllNormal = [1, 0, 0, 0, 0, 0, 0]

#         #default is allnormal
#         PatientTypeProbabilities = PatientTypeProbabilitiesAllNormal
#         patientTypeDistributionOrigin = PatientTypeDistributionOrigin.ALLNORMAL

#         #/ <summary>
        

#         #/ Specifies the probability per patient type of the different lead times [0,14]
        


#         #/ <summary>
#         #/ Check whether a unit is compatible with a request based on the PatientTypeCostMatrix above
#         #/ </summary>
#         #/ <param name="unit"></param>
#         #/ <param name="request"></param>
#         #/ <returns></returns>
#         @staticmethod
#         def IsCompatible(unit, request):
#             type_index = int(request.PatientType)
#             k = 0
#             while k < self.Blood.Antigens.Length:
#                 if RBCallocation.Request.PatientTypeCostMatrix[type_index][k] == -1 and unit.Blood.Vector[k] == 1 and request.Blood.Vector[k] == 0:
#                     return 0
#                 k += 1
#             return 1


#         #/ <summary>
#         #/ Patient's request for blood
#         #/ </summary>
#         #/ <param name="blood"></param>
#         #/ <param name="day_needed"></param>
#         #/ <param name="day_available"></param>
#         #/ <param name="num_units"></param>
#         #/ <param name="patient_type"></param>
#         def __init__(self, blood, day_needed, day_available, num_units, patient_type):
#             #instance fields found by C# to Python Converter:
#             self.DayNeeded = 0
#             self.DayAvailable = 0
#             self.num_units = 0
#             self.Reciprocal = 0
#             self.PatientType = 0
#             self.Blood = None
#             self.MaximumMismatchPenalty = 100
#             self.HospitalIndex = 0
#             self.AllocatedAtDist = False
#             self.Index = 0

#             self.Blood = blood
#             self.DayAvailable = day_available
#             self.DayNeeded = day_needed
#             self.num_units = num_units
#             self.Reciprocal = 1 / float(num_units)
#             self.PatientType = patient_type


#         @staticmethod
#         def FromCsv(csvLine):
#             splitted = csvLine.split(',').ToList()
#             day_needed = int.Parse(splitted[0])
#             day_available = int.Parse(splitted[1])
#             num_units = int.Parse(splitted[2])
#             patient_type = None
#             if splitted[3] == "Wu45":
#                 patient_type = PatientType.WU45
#             elif splitted[3] != "Other":
#                 patient_type = Enum.Parse(typeof(PatientType), splitted[3])
#             else:
#                 patient_type = Enum.Parse(typeof(PatientType), splitted[3])
#             ethnicity = Enum.Parse(typeof(Ethnicity), splitted[4])

#             vector = Array.ConvertAll(splitted.GetRange(5, 17).ToArray(), int.Parse)
#             blood = Blood(vector, ethnicity)
#             return Request(blood, day_needed, day_available, num_units, patient_type)

#         @staticmethod
#         def ToCSVString(r):
#             result = [None for _ in range(5 + 17)]
#             result[0] = str(r.DayNeeded)
#             result[1] = str(r.DayAvailable)
#             result[2] = str(r.NumUnits)
#             result[3] = str(r.PatientType)
#             result[4] = str(r.Blood.Ethnicity)
#             for i in range(0, 17):
#                 result[5 + i] = str(r.Blood.Vector[i])
#             return string.Join(",", result)


#         #/ <summary>
#         #/ Sample a random patient request
#         #/ </summary>
#         #/ <param name="Program.Random"></param>
#         #/ <param name="day_needed"></param>
#         #/ <returns></returns>
#         @staticmethod
        
# }



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
#             self.PatientEthnicityDist[1] = self.Requests.Where(lambda r : r.Blood.Ethnicity == Ethnicity.African).Count() / float(totalNumberOfPatients)
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