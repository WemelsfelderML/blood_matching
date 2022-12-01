from gurobipy import *
import numpy as np
import pandas as pd
import pickle
# import csv
import time
import math
# import os

from blood import *

class Simulation():
    
    def __init__(self, SETTINGS, PARAMS, num_bloodgroups, max_request_size, e):

        # this.InitLength = init_length;
        # this.InventorySize = inventory_size;
        # this.Duration = duration;
        # this.inventory = new Inventory(this.InventorySize);

        # this.HandledRequests = new HandledRequest[demand.GetTotalNumberOfRequests(duration, init_length)];
        
        self.num_bloodgroups = num_bloodgroups
        self.max_request_size = max_request_size
        # self.state_size =  + (num_bloodgroups * max_request_size)

        self.supply_scenario = pd.read_csv(SETTINGS.home_dir + f"supply/{SETTINGS.supply_size}/cau{round(SETTINGS.donor_eth_distr[0]*100)}_afr{round(SETTINGS.donor_eth_distr[1]*100)}_asi{round(SETTINGS.donor_eth_distr[2]*100)}_{e}.csv")
        self.demand_scenario = pd.read_csv(SETTINGS.home_dir + f"demand/{SETTINGS.avg_daily_demand}/{SETTINGS.test_days + SETTINGS.init_days}/{SETTINGS.demand_scenario}_{e}.csv")

        self.inventory = (num_bloodgroups * PARAMS.max_age)

        self.supply_index = 0

        if SETTINGS.strategy == "relimm":
            self.mismatch_weights = PARAMS.relimm_weights
        elif SETTINGS.strategy == "patgroups":
            self.mismatch_weights = PARAMS.patgroup_weights


    def fill_initial_inventory(self, SETTINGS, PARAMS):
        
        inventory = []

        # TODO: now we sample inventory for each day with uniform distribution, change this to sampling more young units and less old units.
        n_products = round(SETTINGS.inventory_size / PARAMS.max_age)

        # Sample supply for each age upto maximum age.
        for age in range(PARAMS.max_age):
            inventory += self.sample_supply_single_day(PARAMS, n_products, age)

        return inventory


    def update_inventory(self, SETTINGS, PARAMS, inventory, x):

        xi = x.sum(axis=1)
        remove = []
        for i in range(len(inventory)):
            # Remove issued products from inventory.
            if xi[i] >= 1:
                remove.append(i)
            # Remove all products that will be outdated at the end of this day.
            elif inventory[i].age >= 34:
                remove.append(i)
            # Increase the age of all remaining products.
            else:
                inventory[i].age += 1

        inventory = [i for i in inventory if i not in remove]
        
        # Generate new supply.
        inventory += self.sample_supply_single_day(PARAMS, max(0, SETTINGS.inventory_size - len(inventory)))       

        return inventory


    def sample_supply_single_day(self, PARAMS, n_products, age = 0):

        # Select the next part of the supply scenario.
        data = self.supply_scenario.iloc[self.supply_index : self.supply_index + n_products]
        self.supply_index += n_products

        inventory = []
        for i in data.index:
            inventory.append(Blood(PARAMS, ethnicity = data.loc[i,"Ethnicity"], major = vector_to_major([data.loc[i,a] for a in PARAMS.major]), minor = [data.loc[i,a] for a in PARAMS.minor], age = age))

        return inventory


    def sample_requests_single_day(self, PARAMS, day = 0):

        # Select the part of the demand scenario belonging to the given day.
        # TODO: checken of het klopt met "Day Available" en "Day Needed"
        data = self.demand_scenario.loc[self.demand_scenario["Day Available"] == day]

        requests = []
        for i in data.index:
            requests.append(Blood(PARAMS, ethnicity = data.loc[i,"Ethnicity"], patgroup = data.loc[i,"Patient Type"], major = vector_to_major([data.loc[i,a] for a in PARAMS.major]), minor = [data.loc[i,a] for a in PARAMS.minor], num_units = data.loc[i,"Num Units"], issuing_day = data.loc[i,"Day Needed"]))

        return requests


# TODO
def minrar(SETTINGS, PARAMS):

    num_bloodgroups = 2**len(PARAMS.major + PARAMS.minor)
    max_request_size = len(PARAMS.request_num_units_probabilities["Other"])

    for e in range(SETTINGS.episodes):
        print(f"\nEpisode: {e}")

        sim = Simulation(SETTINGS, PARAMS, num_bloodgroups, max_request_size, e)
        inventory = sim.fill_initial_inventory(SETTINGS, PARAMS)
        requests = []

        df = SETTINGS.create_output_file(PARAMS, e)

        for day in range(SETTINGS.init_days + SETTINGS.test_days):
            print(f"\nDay {day}")

            requests = [r for r in requests if r.issuing_day >= day]
            requests += sim.sample_requests_single_day(PARAMS, day=day)

            if SETTINGS.strategy == "relimm":
                model, calc_time = allocate_minrar_relimm(SETTINGS, PARAMS, inventory, requests, day)
            elif SETTINGS.strategy == "patgroups":
                model, calc_time = allocate_minrar_patgroups(SETTINGS, PARAMS, inventory, requests, day)

            df.loc[day,"gurobi status"] = model.status
            df.loc[day,"calc time"] = calc_time
            if model.status == 2:   # If an optimal solution was found.
                df, x, y, z = extract_model_solution(SETTINGS, PARAMS, df, model, inventory, requests, e, day)
                if day >= SETTINGS.init_days:
                    df = log_results(SETTINGS, PARAMS, df, model, inventory, requests, day, x=x, y=y)
                inventory = sim.update_inventory(SETTINGS, PARAMS, inventory, x)
                
            elif day >= SETTINGS.init_days:
                print(PARAMS.status_code[model.status])
                df = log_results(SETTINGS, PARAMS, df, model, inventory, requests, day)

            else:
                print(PARAMS.status_code[model.status])

        df.to_csv(SETTINGS.generate_filename("results", "ilp") + f"_{e}.csv", sep=',', index=True)
        

def allocate_minrar_relimm(SETTINGS, PARAMS, inventory, requests, day):

    start = time.perf_counter()

    ################
    ## PARAMETERS ##
    ################

    model = Model(name="model")
    if SETTINGS.show_gurobi_output == False:
        model.Params.LogToConsole = 0
    model.setParam('Threads', SETTINGS.gurobi_threads)
    model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

    I = {i : inventory[i] for i in range(len(inventory))}               # Set of all inventory products.
    R = {r : requests[r] for r in range(len(requests))}                 # Set of all requests.

    antigens = PARAMS.major + PARAMS.minor
    A = {antigens[k] : k for k in range(len(antigens))}   # Set of all antigen indices.
    A_minor = {k : A[k] for k in PARAMS.minor}
    A_no_Fyb = {k : A[k] for k in antigens if k != "Fyb"}

    # TODO extend to patgroups scenario where matches can also be incompatible on minor antigens
    # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
    C = compute_compatibility(PARAMS, I, R)

    # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
    t = [1 - min(1, requests[r].issuing_day - day) for r in R.keys()]

    w = PARAMS.relimm_weights   # Relative immunogenicity weights for each antigen.
    max_age = PARAMS.max_age    # Maximum shelf life of products in inventory    


    ###############
    ## VARIABLES ##
    ###############

    # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
    # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
    # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.

    x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
    y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
    z = model.addVars(len(R), len(A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

    model.update()


    #################
    ## CONSTRAINTS ##
    #################


    # Force y[r] to 1 if not all requested units are satisfied.
    model.addConstrs(R[r].num_units * y[r] + quicksum(x[i,r] for i in I.keys()) >= R[r].num_units for r in R.keys())

    # For each inventory product i, ensure that i can not be issued more than once.
    model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())

    # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on the major antigens.
    model.addConstrs(x[i,r] <= C[i,r] for i in I.keys() for r in R.keys())

    # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
    model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in A_no_Fyb.values())

    # A request can only be mismatched on Fyb if it is positive for Fya.
    model.addConstrs(quicksum(x[i,r] * I[i].vector[A["Fyb"]] * (1 - R[r].vector[A["Fyb"]]) * R[r].vector[A["Fya"]] for i in I.keys()) <= z[r,A["Fyb"]] * R[r].num_units for r in R.keys())    

    # quicksum(z[r,k] * w[k] for r in R.keys() for k in A)                                                                    # Mismatches on minor antigens.
    # quicksum(x[i,r] * w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for i in I.keys() for r in R.keys() for k in A_minor)    # Minor antigen substitution.
    # quicksum(- 1 * math.exp(-4.852 * (max_age - I[i].age - 1) / max_age) * x[i,r] for i in I.keys() for r in R.keys())      # FIFO penalties.
    # quicksum((I[i].get_usability(PARAMS) - R[r].get_usability(PARAMS)) * x[i,r] for i in I.keys() for r in R.keys())        # Product usability on major antigens.


    ################
    ## OBJECTIVES ##
    ################

    # Assign a higher shortage penalty to requests with today as their issuing date.
    model.setObjectiveN(expr = quicksum(y[r] * ((3 * t[r]) + 1) for r in R.keys()), index=0, priority=1, name="shortages") 
    model.setObjectiveN(expr = quicksum(
                                    quicksum(z[r,k] * w[k] for k in A.values()) +                                                   # Mismatches on minor antigens.
                                    quicksum(
                                        quicksum(x[i,r] * w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for k in A_minor.values())   # Minor antigen substitution.
                                        + (- 1 * math.exp(-4.852 * (max_age - I[i].age - 1) / max_age) * x[i,r])                    # FIFO penalties.
                                        + ((I[i].get_usability(PARAMS) - R[r].get_usability(PARAMS)) * x[i,r])                      # Product usability on major antigens.
                                    for i in I.keys()) 
                                for r in R.keys()), index=1, priority=0, name="other")

    # Minimize the objective functions.
    model.ModelSense = GRB.MINIMIZE

    stop = time.perf_counter()
    print(f"model initialization: {(stop - start):0.4f} seconds")


    ##############
    ## OPTIMIZE ##
    ##############

    start = time.perf_counter()
    model.optimize()
    stop = time.perf_counter()
    calc_time = stop - start
    print(f"optimize: {calc_time:0.4f} seconds")

    return model, calc_time


def extract_model_solution(SETTINGS, PARAMS, df, model, inventory, requests, episode, day):

    start = time.perf_counter()

    # Get the values of the model variable as found for the optimal solution.
    x = np.zeros([len(inventory), len(requests)])
    y = np.zeros([len(requests)])
    z = np.zeros([len(requests), len(PARAMS.major + PARAMS.minor)])
    for var in model.getVars():
        name = re.split(r'\W+', var.varName)[0]
        if name == "x":
            index0 = int(re.split(r'\W+', var.varName)[1])
            index1 = int(re.split(r'\W+', var.varName)[2])
            x[index0, index1] = var.X
        if name == "y":
            index0 = int(re.split(r'\W+', var.varName)[1])
            y[index0] = var.X
        if name == "z":
            index0 = int(re.split(r'\W+', var.varName)[1])
            index1 = int(re.split(r'\W+', var.varName)[2])
            z[index0, index1] = var.X

    # Write the values found to local files.
    with open(SETTINGS.generate_filename("results", f"x_{episode}-{day}")+".pickle", "wb") as f:
        pickle.dump(x, f)
    with open(SETTINGS.generate_filename("results", f"y_{episode}-{day}")+".pickle", "wb") as f:
        pickle.dump(y, f)
    with open(SETTINGS.generate_filename("results", f"z_{episode}-{day}")+".pickle", "wb") as f:
        pickle.dump(z, f)
    
    df.loc[day,"nvars"] = sum([np.product(var.shape) for var in [x, y, z]])

    stop = time.perf_counter()
    print(f"loading model results: {(stop - start):0.4f} seconds")

    return df, x, y, z


def log_results(SETTINGS, PARAMS, df, model, inventory, requests, day, x=[], y=[]):

    antigens = PARAMS.major + PARAMS.minor
    ABOD_names = PARAMS.ABOD
    patgroups = PARAMS.patgroups
    ethnicities = ["Caucasian", "African", "Asian"]
    requests = [r for r in requests if r.issuing_day == day]

    df.loc[day,"num patients"] = len(requests)
    df.loc[day,"num units requested"] = sum([r.num_units for r in requests])
    for eth in ethnicities:
        df.loc[day,f"num {eth} patients"] = sum([1 for r in requests if r.ethnicity == eth])
    for p in patgroups:
        df.loc[day,f"num {p} patients"] = sum([1 for r in requests if r.patgroup == p])
        df.loc[day,f"num units requested {p}"] = sum([r.num_units for r in requests if r.patgroup == p])
    
    for i in range(1,5):
        df.loc[day,f"num requests {i} units"] = sum([1 for r in requests if r.num_units == i])

    df.loc[day,"num supplied products"] = sum([1 for i in inventory if i.age == 0])
    for major in ABOD_names:
        df.loc[day,f"num supplied {major}"] = sum([1 for i in inventory if i.major == major and i.age == 0])
        df.loc[day,f"num requests {major}"] = sum([1 for r in requests if r.major == major])
        df.loc[day,f"num {major} in inventory"] = sum([1 for i in inventory if i.major == major])

    # If an optimal solution was found.
    if model.status == 2:

        start = time.perf_counter()

        xi = x.sum(axis=1)  # For each inventory product i∈I, xi[i] = 1 if the product is issued, 0 otherwise.
        xr = x.sum(axis=0)  # For each request r∈R, xr[r] = the number of products issued to this request.
        age_sum = 0
        for r in [r for r in range(len(requests)) if xr[r] > 0]:
            # Get all products from inventory that were issued to request r.
            issued = np.where(x[:,r]==1)[0]
            r = requests[r]

            mismatch = {k:0 for k in antigens}
            for i in [inventory[i] for i in issued]:
                age_sum += i.age
                df.loc[day,f"{i.major} to {r.major}"] += 1
                df.loc[day,f"{i.ethnicity} to {r.ethnicity}"] += 1
            
                # Get all antigens k where product i and request r have a mismatch.
                for k in [antigens[k] for k in range(len(antigens)) if i.vector[k] > r.vector[k]]:
                    # Fy(a-b-) should only be matched on Fy(a), not on Fy(b). -> Fy(b-) only mismatch when Fy(a+)
                    if (k != "Fyb") or (r.vector[antigens.index("Fya")] == 1):    
                        mismatch[k] = 1
                        df.loc[day,[f"num mismatched units {r.patgroup} {k}"]] += 1

            for k in antigens:
                df.loc[day,[f"num mismatches {r.patgroup} {k}"]] += mismatch[k]
                df.loc[day,[f"num mismatches {r.ethnicity} {k}"]] += mismatch[k]

        df.loc[day,f"avg issuing age"] = age_sum / sum(xi)

        for i in [inventory[i] for i in range(len(inventory)) if (xi[i] == 0) and (inventory[i].age >= 34)]:
            df.loc[day,"num outdates"] += 1
            df.loc[day,f"num outdates {i.major}"] += 1

        for r in [requests[r] for r in range(len(requests)) if y[r] == 1]:
            df.loc[day,"num shortages"] += 1
            df.loc[day,f"num shortages {r.major}"] += 1
            df.loc[day,f"num shortages {r.patgroup}"] += 1
            df.loc[day,f"num {r.patgroup} {r.num_units - xr[r]} units short"] += 1
        
        stop = time.perf_counter()
        print(f"writing results to dataframe: {(stop - start):0.4f} seconds")

    
    return df