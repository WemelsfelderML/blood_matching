import pandas as pd

from blood import *

class Hospital():
    
    def __init__(self, SETTINGS, PARAMS, htype, e):

        self.name = f"{htype[:3]}_{e}"

        self.htype = htype
        self.avg_daily_demand = SETTINGS.avg_daily_demand[htype]
        self.inventory_size = SETTINGS.inv_size_factor * self.avg_daily_demand

        self.demand_data = pd.read_csv(SETTINGS.home_dir + f"demand/{self.avg_daily_demand}/{SETTINGS.test_days + SETTINGS.init_days}/{htype}_{e}.csv")

        # Sample supply for each age upto maximum age.
        # TODO maybe have inventory and requests already be a index-product dictionary as used in MINRAR?
        self.inventory = []
        self.requests = []

    def update_inventory(self, SETTINGS, PARAMS, x, day):

        I = {i : self.inventory[i] for i in range(len(self.inventory))}
        remove = []

        # Remove all products form inventory that were issued to requests with today as their issuing date
        for r in [r for r in range(len(self.requests)) if self.requests[r].issuing_day == day]:
            remove += list(np.where(x[:,r]==1)[0])

        # If a product will be outdated at the end of this day, remove from inventory, otherwise increase its age.
        for i in I.keys():
            if I[i].age >= (PARAMS.max_age-1):
                remove.append(i)
            else:
                I[i].age += 1
 
        self.inventory = [I[i] for i in I.keys() if i not in remove]

        # Return the number of products to be supplied, in order to fill the inventory upto its maximum capacity.
        return max(0, self.inventory_size - len(self.inventory))


    def sample_requests_single_day(self, PARAMS, day = 0):

        # Select the part of the demand scenario belonging to the given day.
        # TODO: checken of het klopt met "Day Available" en "Day Needed"
        data = self.demand_data.loc[self.demand_data["Day Available"] == day]

        requests = []
        for i in data.index:
            requests.append(Blood(PARAMS, ethnicity = data.loc[i,"Ethnicity"], patgroup = data.loc[i,"Patient Type"], major = vector_to_major([data.loc[i,a] for a in PARAMS.major]), minor = [data.loc[i,a] for a in PARAMS.minor], num_units = data.loc[i,"Num Units"], issuing_day = data.loc[i,"Day Needed"]))

        self.requests += requests