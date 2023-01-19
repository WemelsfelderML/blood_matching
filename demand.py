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



class Demand:

    
    def __init__(self, weekday_index, avg_daily_demand):

        national_demand = [[2.160824e+03, 7.144959e-02], [1.940279e+03, 8.390841e-02], [1.980793e+03, 7.903061e-02], [1.954597e+03, 7.880648e-02], [2.029441e+03, 7.298936e-02], [319.5225256, 0.2143665], [308.0818815, 0.2747322]]
        avg_national_demand = (national_demand[0][0] + national_demand[1][0] + national_demand[2][0] + national_demand[3][0] + national_demand[4][0] + national_demand[5][0] + national_demand[6][0]) / 7

        # sizefactor = number of identical subdistributions which together form exactly the national_demand distribution
        sizefactor = avg_national_demand / float(avg_daily_demand)

        national_mean = national_demand[weekday_index][0]
        national_cv = national_demand[weekday_index][1]
        national_stdev = national_cv * national_mean

        # compute mean and stdev for subdistribution (= scaled down distribution)
        mean = national_mean / sizefactor
        stdev = national_stdev / float(math.sqrt(sizefactor))

        # procedure below is based on Source: https://www.win.tue.nl/~iadan/alqt/fit.pdf
        # compute parameters for fitting
        self._a = (stdev / mean) ** 2 - 1 / mean

        # mixed geometric distribution
        if self._a >= 1:
            self._p1 = (mean * (1 + self._a + math.sqrt(self._a * self._a - 1))) / (2 + mean * (1 + self._a + math.sqrt(self._a * self._a - 1)))
            self._p2 = (mean * (1 + self._a - math.sqrt(self._a * self._a - 1))) / (2 + mean * (1 + self._a - math.sqrt(self._a * self._a - 1)))
            self._q1 = 1 / (1 + self._a + math.sqrt(self._a * self._a - 1))
            self._q2 = 1 / (1 + self._a - math.sqrt(self._a * self._a - 1))

        # mixed negative binomial distribution
        else:
            self._k = math.floor(1 / self._a)
            self._q = ((self._k + 1) * self._a - math.sqrt((self._k + 1) * (1 - self._a * self._k))) / (1 + self._a)
            self._p = mean / ((self._k + 1) - self._q + mean)

    # Sample a number of units demanded based on parameters as computed in the constructor of this class
    # Source: https://www.win.tue.nl/~iadan/alqt/fit.pdf
    def sample_number_of_units(self):
        #For sampling the geometric distribution we use (1-p) as parameter and subtract one to account for the difference in 
        #definitions for the geometric distribution. See https://en.wikipedia.org/wiki/Geometric_distribution#:~:text=The%20geometric%20distribution%20gives%20the,%2C%203%2C%20
        #Paper uses the definition on the right, whereas the implementation uses the definition on the left.

        if self._a >= 1: #sample geometric
            if random.random() < self._q1:
                return self.sample_geometric(1 - self._p1) - 1
            else:
                return self.sample_geometric(1 - self._p2) - 1
        else: #sample negative binomial (repeated geometric sampling)
            sum = 0
            if random.random() < self._q:
                i = 0
                while i < self._k:
                    sum += self.sample_geometric(1 - self._p) - 1
                    i += 1
            else:
                i = 0
                while i < self._k+1:
                    sum += self.sample_geometric(1 - self._p) - 1
                    i += 1
            return sum

    # Generate a list of random requests according to the given distribution
    def sample_request_for_day(self, SETTINGS, PARAMS, day_index, df, htype):
        
        # Keep sampling new requests for today until the required number of requests is reached.
        num_units_requested = self.sample_number_of_units()
        num_units = 0
        while num_units < num_units_requested:
            r = self.get_random_request(SETTINGS, PARAMS, day_index, htype)
            df.loc[len(df)] = [r.day_needed, r.day_available, r.num_units, r.patgroup, r.blood.ethnicity] + r.blood.vector
            num_units += r.num_units

        return df

    def get_random_request(self, SETTINGS, PARAMS, day_needed, htype):

        # Determine the patient group, lead time, number of units and ethnicity of the patient request.
        patgroup = random.choices(PARAMS.patgroups, weights = [PARAMS.patgroup_distr[htype][pg] for pg in PARAMS.patgroups], k=1)[0]
        lead_time = random.choices(range(14), weights = PARAMS.request_lead_time_probabilities[patgroup], k=1)[0]
        num_units = random.choices(range(1,5), weights = PARAMS.request_num_units_probabilities[patgroup], k=1)[0]

        if patgroup == "SCD":
            ethnicity = "African"
        else:
            ethnicity = "Caucasian"

        # Create a Request instance using the generated information.
        return Request(day_needed, max(0, day_needed - lead_time), num_units, patgroup, Blood(PARAMS, ethnicity))


    # Sample a geometric distribution with parameter p (mean = 1/p)
    # This method is a replacement for MathNet.Numerics.Distributions.Geometric(p)
    def sample_geometric(self, p):
        return int(math.ceil(math.log(1 - random.random()) / math.log(1 - p)))




# Generate multiple demand scenarios with the same parameters.
def generate_demand(SETTINGS, PARAMS, htype, avg_daily_demand):

    duration = SETTINGS.test_days + SETTINGS.init_days

    # Make sure all required folders exist.
    path = SETTINGS.home_dir + "demand"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"demand/{avg_daily_demand}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    # Find already existing demand files of the chosen size and duration, and make sure not to overwrite them.
    i = 0
    while os.path.exists(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{htype}_{i}.csv"):
        i += 1

    # For every episode in the given range, generate requests for all days of the simulation.
    for _ in range(SETTINGS.episodes[0],SETTINGS.episodes[1]):

        print(f"Generating demand '{htype}_{i}'.")

        # initialize distributions for each day of the week
        daily_distributions = []
        for weekday_index in range(7):
            # daily_distributions.append(Demand(SETTINGS, weekday_index, avg_daily_demand))
            daily_distributions.append(Demand(weekday_index, avg_daily_demand))

        df = pd.DataFrame(columns = ["Day Needed", "Day Available", "Num Units", "Patient Type", "Ethnicity"] + PARAMS.major + PARAMS.minor)

        requests = []
        for day_index in range(duration):
            df = daily_distributions[day_index%7].sample_request_for_day(SETTINGS, PARAMS, day_index, df, htype)

        df.to_csv(SETTINGS.home_dir + f"demand/{avg_daily_demand}/{duration}/{htype}_{i}.csv", index=False)

        i += 1



