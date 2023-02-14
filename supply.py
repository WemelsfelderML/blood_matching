import math
import os
import pandas as pd

from blood import *


# Generate a given number of supply files, where each file contains enough supply for one simulation episode.
def generate_supply(SETTINGS, PARAMS):

    size = SETTINGS.supply_size                                                     # Total number of products to be generated.
    eth = SETTINGS.donor_eth_distr                                                  # Proportion of Caucasian, African and Asian donors.
    name = f"cau{round(eth[0]*100)}_afr{round(eth[1]*100)}_asi{round(eth[2]*100)}"  # Name to give the supply file.

    # Make sure all required folders exist.
    path = SETTINGS.home_dir + "supply"
    if os.path.exists(path) == False:
        os.mkdir(path)

    path = SETTINGS.home_dir + f"supply/{size}"
    if os.path.exists(path) == False:
        os.mkdir(path)

    # Find already existing supply files of the chosen size and duration, and make sure not to overwrite them.
    i = 0
    while os.path.exists(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv"):
        i += 1

    # For every episode in the given range, generate enough supply for each simulation.
    for _ in range(SETTINGS.episodes[0],SETTINGS.episodes[1]):

        print(f"Generating supply '{name}_{i}'.")

        # Generate the required number of products and write to a pandas dataframe.
        products = generate_products(SETTINGS, PARAMS, size)
        df = pd.DataFrame(products, columns = PARAMS.major + PARAMS.minor + ["Ethnicity"])
        
        # Shuffle the supplied products.
        index = list(df.index)
        random.shuffle(index)
        df = df.loc[index]

        # # Uncomment the code below to let the initial inventory consist of products of only one major blood type.
        # if sum(SETTINGS.n_hospitals.values()) > 1:
        #     inventory_size = SETTINGS.inv_size_factor_dc * sum([SETTINGS.n_hospitals[htype] * SETTINGS.avg_daily_demand[htype] for htype in SETTINGS.n_hospitals.keys()])
        # else:
        #     inventory_size = SETTINGS.inv_size_factor_hosp * sum([SETTINGS.n_hospitals[htype] * SETTINGS.avg_daily_demand[htype] for htype in SETTINGS.n_hospitals.keys()])
        # products = []
        # for _ in range(inventory_size):
        #     ip = Blood(PARAMS, ethnicity="Caucasian", major="AB+")    # Change the 'major' argument to the major blood group to use for the initial inventory.
        #     products.append(ip.vector + [ip.ethnicity])
        # df = pd.concat([pd.DataFrame(products, columns = PARAMS.major+PARAMS.minor+["Ethnicity"]), df.iloc[inventory_size:]])

        df["Index"] = range(len(df))
        
        # Write the dataframe to a csv file.
        df.to_csv(SETTINGS.home_dir + f"supply/{size}/{name}_{i}.csv", index=False)

        i += 1

# Generate a list of products with a specific ethnic distribution and a specific ABODistribution, in random order
def generate_products(SETTINGS, PARAMS, size):

    products = []
    majors_sampled = {major : 0 for major in PARAMS.ABOD}

    # Sample African donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[1])):
        ip = Blood(PARAMS, ethnicity="African")
        products.append(ip.vector + [ip.ethnicity])
        majors_sampled[ip.major] += 1

    # Sample Asian donors.
    for _ in range(round(size * SETTINGS.donor_eth_distr[2])):
        ip = Blood(PARAMS, ethnicity="Asian")
        products.append(ip.vector + [ip.ethnicity])
        majors_sampled[ip.major] += 1

    # for _ in range(round(size * SETTINGS.donor_eth_distr[0])):
    #     ip = Blood(PARAMS, "Caucasian")
    #     products.append(ip.vector + [ip.ethnicity])

    # For each major blood group determine how many products should be additionally sampled, to make sure that the overall list has the correct ABOD distribution
    for major in PARAMS.ABOD:
        num_to_sample = round(PARAMS.donor_ABOD_distr[major] * size) - majors_sampled[major]
        
        for _ in range(max(0,num_to_sample)):
            ip = Blood(PARAMS, ethnicity="Caucasian", major=major)
            products.append(ip.vector + [ip.ethnicity])

    return products
