import math
import os
import pandas as pd

from blood import *

class Unit:

    def __init__(self, blood):

        self.blood = blood
        self.age = 0
        self.shelf_life = 35

    # def get_fifo_cost(self):
    #     return float((1 - math.exp(-4.852 * (self.shelf_life - self.Age - 1) / self.shelf_life)))

    # def GetUsability(self, antigen_set):
    #     return self.Blood.Usability(antigen_set)

    # def GetRelativeOpportunityLoss(self, request, antigen_set):
    #     unit_usab = self.GetUsability(antigen_set)
    #     request_usab = request.Blood.Usability(antigen_set)
    #     rol = (unit_usab - request_usab) / unit_usab
    #     return float((1 - math.exp(-4.852 * rol)))

    # def FromCsv(csvLine):
    #     splitted = csvLine.split(',')

    #     #parse the presence of all 17 antigens
    #     vector = Array.ConvertAll(splitted.Take(17).ToArray(), int.Parse)

    #     #if a ethnicity is given, then also parse this. Otherwise use Caucasian as default
    #     eth = Ethnicity.Caucasian
    #     if len(splitted)>17:
    #         eth = Enum.Parse(typeof(Ethnicity), splitted[17])


    #     blood = Blood(vector, eth)
    #     return Unit(blood, 0)

    # def ToCsv(u):
    #     return string.Join(",", u.Blood.Vector.Select(lambda x : str(x)).ToArray())+","+str(u.Blood.Ethnicity)


# Generate multiple scenarios with the same parameters
def generate_supply(SETTINGS):

    size = SETTINGS.supply_size
    name = SETTINGS.name

    path = SETTINGS.home_dir + "supply"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"supply/{size}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    i = 0
    while os.path.exists(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv"):
        i += 1

    for _ in range(SETTINGS.episodes):

        print(f"Generating supply '{name}_{i}'.")

        units = generate_units(SETTINGS, size)
        df = pd.DataFrame(units, columns = SETTINGS.major + SETTINGS.minor + ["Ethnicity"])
        
        index = list(df.index)
        random.shuffle(index)
        df = df.loc[index]
        
        df.to_csv(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv", index=False)

        i += 1

# Generate a list of units with a specific ethnic distribution and a specific ABODistribution, in random order
def generate_units(SETTINGS, size):

    units = []

    majors_sampled = {major : 0 for major in SETTINGS.ABOD}

    # Sample African donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[1])):
        unit = Unit(Blood(SETTINGS, "African"))
        units.append(unit.blood.vector + [unit.blood.ethnicity])
        majors_sampled[unit.blood.major] += 1

    # Sample Asian donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[2])):
        unit = Unit(Blood(SETTINGS, "Asian"))
        units.append(unit.blood.vector + [unit.blood.ethnicity])
        majors_sampled[unit.blood.major] += 1

    # For each major blood group determine how many units should be additionally sampled to make sure that the overall list has the correct ABOD distribution
    for major in SETTINGS.ABOD:
        num_to_sample = round(SETTINGS.donor_ABOD_distr[major] * size) - majors_sampled[major]
        
        for _ in range(max(0,num_to_sample)):
            unit = Unit(Blood(SETTINGS, "Caucasian", major=major))
            units.append(unit.blood.vector + [unit.blood.ethnicity])

    return units


