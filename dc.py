import pandas as pd

from blood import *

class Distribution_center():
    
    def __init__(self, SETTINGS, PARAMS, hospitals, e):

        self.name = f"dc_{e}"

        self.supply_data = pd.read_csv(SETTINGS.home_dir + f"supply/{SETTINGS.supply_size}/cau{round(SETTINGS.donor_eth_distr[0]*100)}_afr{round(SETTINGS.donor_eth_distr[1]*100)}_asi{round(SETTINGS.donor_eth_distr[2]*100)}_{e}.csv")
        self.supply_index = 0

        if len(hospitals) > 1:

            self.inventory_size = SETTINGS.inv_size_factor * sum([hospital.avg_daily_demand for hospital in hospitals])

            # Sample supply for each age upto maximum age.
            inventory = []
            n_products = round(self.inventory_size / PARAMS.max_age)
            for age in range(PARAMS.max_age):
                inventory += self.sample_supply_single_day(PARAMS, n_products, age)

            self.inventory = inventory


    def update_inventory(self, SETTINGS, PARAMS, x, day):

        I = {i : self.inventory[i] for i in range(len(self.inventory))}
        remove = []

        # Remove all products form inventory that were issued to requests with today as their issuing date
        xi = x.sum(axis=1)  # For each inventory product i∈I, xi[i] = 1 if the product is transported, 0 otherwise.
        remove += [i for i in I.keys() if xi[i] >= 1]

        # If a product will be outdated at the end of this day, remove from inventory, otherwise increase its age.
        for i in I.keys():
            if I[i].age >= (PARAMS.max_age-1):
                remove.append(i)
            else:
                I[i].age += 1

        self.inventory = [I[i] for i in I.keys() if i not in remove]
        
        # Generate new supply.
        self.inventory += self.sample_supply_single_day(PARAMS, max(0, self.inventory_size - len(self.inventory)))


    def sample_supply_single_day(self, PARAMS, n_products, age = 0):

        # Select the next part of the supply scenario.
        data = self.supply_data.iloc[self.supply_index : self.supply_index + n_products]
        self.supply_index += n_products

        supply = []
        for i in data.index:
            supply.append(Blood(PARAMS, ethnicity = data.loc[i,"Ethnicity"], major = vector_to_major([data.loc[i,a] for a in PARAMS.major]), minor = [data.loc[i,a] for a in PARAMS.minor], age = age))

        return supply
