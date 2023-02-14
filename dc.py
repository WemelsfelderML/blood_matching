import pandas as pd

from blood import *

class Distribution_center():
    
    def __init__(self, SETTINGS, PARAMS, hospitals, e):

        # Name for the distribution center (currently not used)
        self.name = f"dc_{e}"

        # Read the supply that was generated using SETTINGS.mode = "supply"
        self.supply_data = pd.read_csv(SETTINGS.home_dir + f"supply/{SETTINGS.supply_size}/cau{round(SETTINGS.donor_eth_distr[0]*100)}_afr{round(SETTINGS.donor_eth_distr[1]*100)}_asi{round(SETTINGS.donor_eth_distr[2]*100)}_{e}.csv")
        
        # Keep track of the supply index to know which item of the supply data to read next.
        self.supply_index = 0

        # In the multi-hospital scenario, the distribution center also has its own inventory.
        if len(hospitals) > 1:

            self.inventory_size = SETTINGS.inv_size_factor_dc * sum([hospital.avg_daily_demand for hospital in hospitals])

            # Initialize the inventory with products from the supply data, where the product's age is uniformly distributed between 0 and the maximum shelf life.
            inventory = []
            n_products = round(self.inventory_size / PARAMS.max_age)
            for age in range(PARAMS.max_age):
                inventory += self.sample_supply_single_day(PARAMS, n_products, age)

            self.inventory = inventory


    # Update the distribution centers's inventory at the end of a day in the simulation.
    def update_inventory(self, SETTINGS, PARAMS, x, day):

        I = {i : self.inventory[i] for i in range(len(self.inventory))}
        remove = []

        # Remove all products from inventory that will be shipped to a hospital.
        xi = x.sum(axis=1)
        remove += [i for i in I.keys() if xi[i] >= 1]

        # If a product will be outdated at the end of this day, remove from inventory, otherwise increase its age.
        for i in I.keys():
            if I[i].age >= (PARAMS.max_age-1):
                remove.append(i)
            else:
                I[i].age += 1

        self.inventory = [I[i] for i in I.keys() if i not in remove]
        
        # Supply the inventory upto its capacity with new products from the supply data.
        self.inventory += self.sample_supply_single_day(PARAMS, max(0, self.inventory_size - len(self.inventory)))


    # Read the required number of products from the supply data and add these products to the distribution center's inventory.
    def sample_supply_single_day(self, PARAMS, n_products, age = 0):

        # Select the next part of the supply scenario.
        data = self.supply_data.iloc[self.supply_index : self.supply_index + n_products]
        self.supply_index += n_products

        # Transform the newly received supply, as read from the data file, to instances of the Blood class.
        supply = []
        for i in data.index:
            supply.append(Blood(PARAMS, index = data.loc[i,"Index"], ethnicity = data.loc[i,"Ethnicity"], major = vector_to_major([data.loc[i,a] for a in PARAMS.major]), minor = [data.loc[i,a] for a in PARAMS.minor], age = age))

        return supply
