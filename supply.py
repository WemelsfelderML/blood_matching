import math
import os
import pandas as pd

from blood import *

class Unit:

    def __init__(self, blood):

        self.blood = blood
        self.age = 0
        self.shelf_life = 35

# Generate multiple scenarios with the same parameters
def generate_supply(SETTINGS, PARAMS):

    size = SETTINGS.supply_size
    eth = SETTINGS.donor_eth_distr
    name = f"cau{round(eth[0]*100)}_afr{round(eth[1]*100)}_asi{round(eth[2]*100)}"

    path = SETTINGS.home_dir + "supply"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"supply/{size}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    i = 0
    while os.path.exists(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv"):
        i += 1

    for _ in range(SETTINGS.episodes[0],SETTINGS.episodes[1]):

        print(f"Generating supply '{name}_{i}'.")

        units = generate_units(SETTINGS, PARAMS, size)
        df = pd.DataFrame(units, columns = PARAMS.major + PARAMS.minor + ["Ethnicity"])
        
        index = list(df.index)
        random.shuffle(index)
        df = df.loc[index]

        inventory_size = SETTINGS.inv_size_factor_hosp * sum([SETTINGS.n_hospitals[htype] * SETTINGS.avg_daily_demand[htype] for htype in SETTINGS.n_hospitals.keys()])
        units = []
        for _ in range(inventory_size):
            unit = Unit(Blood(PARAMS, "Caucasian", major="AB+"))
            units.append(unit.blood.vector + [unit.blood.ethnicity])
        df = pd.concat([pd.DataFrame(units, columns = PARAMS.major+PARAMS.minor+["Ethnicity"]), df.iloc[inventory_size:]])
        
        df.to_csv(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv", index=False)

        i += 1

# Generate a list of units with a specific ethnic distribution and a specific ABODistribution, in random order
def generate_units(SETTINGS, PARAMS, size):

    units = []

    majors_sampled = {major : 0 for major in PARAMS.ABOD}

    # Sample African donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[1])):
        unit = Unit(Blood(PARAMS, "African"))
        units.append(unit.blood.vector + [unit.blood.ethnicity])
        majors_sampled[unit.blood.major] += 1

    # Sample Asian donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[2])):
        unit = Unit(Blood(PARAMS, "Asian"))
        units.append(unit.blood.vector + [unit.blood.ethnicity])
        majors_sampled[unit.blood.major] += 1

    # for _ in range(round(size * SETTINGS.donor_eth_distr[0])):
    #     unit = Unit(Blood(PARAMS, "Caucasian"))
    #     units.append(unit.blood.vector + [unit.blood.ethnicity])

    # For each major blood group determine how many units should be additionally sampled to make sure that the overall list has the correct ABOD distribution
    for major in PARAMS.ABOD:
        num_to_sample = round(PARAMS.donor_ABOD_distr[major] * size) - majors_sampled[major]
        
        for _ in range(max(0,num_to_sample)):
            unit = Unit(Blood(PARAMS, "Caucasian", major=major))
            units.append(unit.blood.vector + [unit.blood.ethnicity])

    return units
