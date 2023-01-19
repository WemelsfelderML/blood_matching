import numpy as np
import pandas as pd
import pickle
# import csv
# import os

from blood import *
from hospital import *
from dc import *
from minrar import *

def simulation(SETTINGS, PARAMS):

    optimizer = MINRAR(SETTINGS, PARAMS)

    if sum(SETTINGS.n_hospitals.values()) > 1:

        for e in range(SETTINGS.episodes[0], SETTINGS.episodes[1]):
            print(f"\nEpisode: {e}")

            hospitals = []
            for htype in SETTINGS.n_hospitals.keys():
                hospitals += [Hospital(SETTINGS, PARAMS, htype, (e*SETTINGS.n_hospitals[htype])+i) for i in range(SETTINGS.n_hospitals[htype])]
            dc = Distribution_center(SETTINGS, PARAMS, hospitals, e)

            # Initialize hospital inventories
            for hospital in hospitals:
                n_products = round(hospital.inventory_size / PARAMS.max_age)
                for age in range(PARAMS.max_age):
                    hospital.inventory += dc.sample_supply_single_day(PARAMS, n_products, age)

            df = SETTINGS.initialize_output_dataframe(PARAMS, hospitals, e)
            
            for day in range(SETTINGS.init_days + SETTINGS.test_days):
                print(f"\nDay {day}")

                df = simulate_multiple_hospitals(SETTINGS, PARAMS, optimizer, df, dc, hospitals, e, day)

            df.to_csv(SETTINGS.generate_filename("results") + f"{SETTINGS.strategy}_{'-'.join([str(SETTINGS.n_hospitals[ds]) + ds[:3] for ds in SETTINGS.n_hospitals.keys()])}_{e}.csv", sep=',', index=True)        
            
    else:

        htype = max(SETTINGS.n_hospitals, key = lambda i: SETTINGS.n_hospitals[i])

        for e in range(SETTINGS.episodes[0], SETTINGS.episodes[1]):
            print(f"\nEpisode: {e}")

            hospital = Hospital(SETTINGS, PARAMS, htype, e)
            dc = Distribution_center(SETTINGS, PARAMS, [hospital], e)

            if SETTINGS.line == "on":
                # Initialize hospital inventory
                n_products = round(hospital.inventory_size / PARAMS.max_age)
                for age in range(PARAMS.max_age):
                    hospital.inventory += dc.sample_supply_single_day(PARAMS, n_products, age)

                df = SETTINGS.initialize_output_dataframe(PARAMS, [hospital], e)

                for day in range(SETTINGS.init_days + SETTINGS.test_days):
                    print(f"\nDay {day}")

                    df = simulate_single_hospital(SETTINGS, PARAMS, optimizer, df, dc, hospital, e, day)

                df.to_csv(SETTINGS.generate_filename("results") + f"{SETTINGS.strategy}_{htype[:3]}_{e}.csv", sep=',', index=True)

            elif SETTINGS.line == "off":

                days = range(SETTINGS.init_days + SETTINGS.test_days)

                hospital.inventory = dc.sample_supply_single_day(PARAMS, SETTINGS.supply_size)
                for day in days:
                    hospital.sample_requests_single_day(PARAMS, day=day)

                df = SETTINGS.initialize_output_dataframe(PARAMS, [hospital], e)
                
                model = optimizer.minrar_offline(SETTINGS, PARAMS, hospital, days)
                
                df, x, y, z, a, b = model_output_to_matches(SETTINGS, PARAMS, df, model, [hospital], e)
                for day in days:
                    df = optimizer.log_results(SETTINGS, PARAMS, df, model, hospital, day, x=x, y=y, z=z, a=a, b=b)

                df.to_csv(SETTINGS.generate_filename("results") + f"offline_{SETTINGS.strategy}_{htype[:3]}_{e}.csv", sep=',', index=True)

            else:
                print("Please set the 'line' variable in settings.py to a valid value (either 'off' or 'on').")


def simulate_single_hospital(SETTINGS, PARAMS, optimizer, df, dc, hospital, e, day):

    hospital.requests = [r for r in hospital.requests if r.issuing_day >= day]
    hospital.sample_requests_single_day(PARAMS, day=day)

    model = optimizer.minrar_single_hospital(SETTINGS, PARAMS, hospital, day)

    df.loc[(day, hospital.name),"gurobi status"] = model.status
    # df.loc[(day, hospital.name),"calc time"] = calc_time

    df, x, y, z = model_output_to_matches(SETTINGS, PARAMS, df, model, dc, [hospital], e, day)
    if day >= SETTINGS.init_days:
        df = optimizer.log_results(SETTINGS, PARAMS, df, model, hospital, day, x=x, y=y, z=z)
    supply_size = hospital.update_inventory(SETTINGS, PARAMS, x, day)
    hospital.inventory += dc.sample_supply_single_day(PARAMS, supply_size)

    return df


def simulate_multiple_hospitals(SETTINGS, PARAMS, optimizer, df, dc, hospitals, e, day):

    for hospital in hospitals:
        hospital.requests = [r for r in hospital.requests if r.issuing_day >= day]
        hospital.sample_requests_single_day(PARAMS, day=day)

        for rq in hospital.requests:
            rq.allocated_from_dc = 0

    model = optimizer.minrar_multiple_hospitals(SETTINGS, PARAMS, dc, hospitals, day)
    df, xh, xdc, y, z = model_output_to_matches(SETTINGS, PARAMS, df, model, dc, hospitals, e, day)
 
    allocations_from_dc = np.zeros([len(dc.inventory),len(hospitals)])
    for h in range(len(hospitals)):
        for r in range(len(hospitals[h].requests)):
            if hospitals[h].requests[r].issuing_day == (day + 1):

                issued = np.where(xdc[h][:,r]==1)[0]
                for i in issued:
                    allocations_from_dc[i,h] = 1
                    hospitals[h].requests[r].allocated_from_dc += 1     # total number of products allocated to this request from DC

        if day >= SETTINGS.init_days:
            df = optimizer.log_results(SETTINGS, PARAMS, df, model, hospitals[h], day, x=xh[h], y=y[h], z=z[h])

    # calculate assignment within hospitals knowing allocations from DC
    supply_sizes = [hospitals[h].update_inventory(SETTINGS, PARAMS, xh[h], day) for h in range(len(hospitals))]

    # print("supply sizes", supply_sizes)
    # print("allocations from dc", allocations_from_dc)
    # print(len(dc.inventory))
    # print([len(hospital.inventory) for hospital in hospitals])
    
    # given supply_sizes, allocate products to each of the hospitals to restock them upto their maximum capacity
    model = optimizer.allocate_remaining_supply_from_dc(SETTINGS, PARAMS, day, dc.inventory, hospitals, supply_sizes, allocations_from_dc)
    x = model_output_to_transports(model, len(dc.inventory), len(hospitals))

    for h in range(len(hospitals)):
        hospitals[h].inventory += [dc.inventory[i] for i in range(len(dc.inventory)) if x[i,h] >= 1]

    # Remove all shipped and outdated units from DC inventory, and increase the age of all remaining products.
    dc.update_inventory(SETTINGS, PARAMS, x, day)

    return df


# def simulate_multiple_hospitals(SETTINGS, PARAMS, optimizer, df, dc, hospitals, e, day):
 
#     propagated_requests = []

#     for h in range(len(hospitals)):
#         hospital = hospitals[h]
#         hospital.requests = [r for r in hospital.requests if r.issuing_day >= day]
#         for rq in hospital.requests:
#             rq.allocated_from_dc = 0
#         hospital.sample_requests_single_day(PARAMS, day=day)

#         model = optimizer.minrar_hospital_before_dc_allocation(SETTINGS, PARAMS, day, hospital)
#         df, x, y, z = model_output_to_matches(SETTINGS, PARAMS, df, model, hospital.inventory, hospital.requests, hospital.name, e, day)

#         for r in range(len(hospital.requests)):
#             if y[r] == 0:
#                 hospital.requests[r].best_mismatch_penalty = optimizer.get_mismatch_penalty(SETTINGS, hospital, r, z)
#             rq = hospital.requests[r]
#             if (rq.issuing_day > day) and (rq.best_mismatch_penalty > 0):
#                 if rq.patgroup in ["MDS", "Thal", "AIHA", "ALA", "SCD"]:
#                     propagated_requests.append([rq, h, r])
#                 elif (rq.patgroup == "Wu45") and (y[r] == 1):
#                     propagated_requests.append([rq, h, r])

#     # calculate product allocations from DC and return a list of length 7 (each day of the coming week),
#     # with for each day a set of (hospital, product) pairs to be shipped
#     model = optimizer.allocate_propagated_requests_from_dc(SETTINGS, PARAMS, day, dc.inventory, propagated_requests)
#     df, x, y, z = model_output_to_matches(SETTINGS, PARAMS, df, model, dc.inventory, [pr[0] for pr in propagated_requests], dc.name, e, day)

#     allocations_from_dc = np.zeros([len(dc.inventory),len(hospitals)])
#     for r in range(len(propagated_requests)):
#         if propagated_requests[r][0].issuing_day == (day + 1):

#             issued = np.where(x[:,r]==1)[0]
#             for i in issued:
#                 allocations_from_dc[i,propagated_requests[r][1]] = 1

#         if y[r] == 0:
#             hospitals[propagated_requests[r][1]].requests[propagated_requests[r][2]].allocated_from_dc = 1

#     supply_sizes = []
#     for hospital in hospitals:
#         # calculate assignment within hospitals knowing allocations from DC
#         model = optimizer.mirar_hospital_after_dc_allocation(SETTINGS, PARAMS, day, hospital)
#         df, x, y, z = model_output_to_matches(SETTINGS, PARAMS, df, model, hospital.inventory, hospital.requests, hospital.name, e, day)

#         if day >= SETTINGS.init_days:
#             df = optimizer.log_results(SETTINGS, PARAMS, df, model, hospital, day, x=x, y=y, z=z)
#         supply_sizes.append(hospital.update_inventory(SETTINGS, PARAMS, x, day))

    
#     # given supply_sizes, allocate products to each of the hospitals to restock them upto their maximum capacity
#     model = optimizer.allocate_remaining_supply_from_dc(SETTINGS, PARAMS, day, dc.inventory, hospitals, supply_sizes, allocations_from_dc)

#     x = model_output_to_transports(model, len(dc.inventory), len(hospitals))

#     # I = {i : dc.inventory[i] for i in range(len(dc.inventory))}               # Set of all inventory products.
#     # H = {h : hospitals[h] for h in range(len(hospitals))}               # Set of all inventory products.
#     # bi = [I[i].get_usability(PARAMS, hospitals, antigens=["C", "c", "E", "e", "K", "k", "Fya", "Fyb", "Jka", "Jkb"]) for i in I.keys()]
    
#     for h in range(len(hospitals)):
#         hospitals[h].inventory += [dc.inventory[i] for i in range(len(dc.inventory)) if x[i,h] >= 1]

#     # Remove all shipped and outdated units from DC inventory, and increase the age of all remaining products.
#     dc.update_inventory(SETTINGS, PARAMS, x, day)
#     # TODO: log outdates in dc

#     return df


def model_output_to_matches(SETTINGS, PARAMS, df, model, dc, hospitals, episode, day=0):

    # start = time.perf_counter()

    if sum(SETTINGS.n_hospitals.values()) > 1:

        xh = [np.zeros([len(hospital.inventory), len(hospital.requests)]) for hospital in hospitals]
        xdc = [np.zeros([len(dc.inventory), len(hospital.requests)]) for hospital in hospitals]
        y = [np.zeros([len(hospital.requests)]) for hospital in hospitals]
        z = [np.zeros([len(hospital.requests), len(PARAMS.major + PARAMS.minor)]) for hospital in hospitals]

        for h in range(len(hospitals)):
            for var in model.getVars():
                var_name = re.split(r'\W+', var.varName)[0]
                if var_name == f"xh{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xh[h][index0, index1] = var.X
                if var_name == f"xdc{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xdc[h][index0, index1] = var.X
                if var_name == f"y{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    y[h][index0] = var.X
                if var_name == f"z{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    z[h][index0, index1] = var.X

        nvars = sum([sum([np.product(var.shape) for var in [xh[h], xdc[h], y[h], z[h]]]) for h in range(len(hospitals))])
        print("nvars:",nvars)
        df.loc[(day,hospitals[h].name),"nvars"] = nvars

        return df, xh, xdc, y, z

    else:

        # Get the values of the model variable as found for the optimal solution.
        x = np.zeros([len(hospitals[0].inventory), len(hospitals[0].requests)])
        y = np.zeros([len(hospitals[0].requests)])
        z = np.zeros([len(hospitals[0].requests), len(PARAMS.major + PARAMS.minor)])

        for var in model.getVars():
            var_name = re.split(r'\W+', var.varName)[0]
            if var_name == "x":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                x[index0, index1] = var.X
            if var_name == "y":
                index0 = int(re.split(r'\W+', var.varName)[1])
                y[index0] = var.X
            if var_name == "z":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                z[index0, index1] = var.X

        if SETTINGS.line == "off":
            # o = np.zeros([len(inventory)])
            a = np.zeros([len(hospitals[0].inventory), SETTINGS.init_days + SETTINGS.test_days])
            b = np.zeros([len(hospitals[0].inventory), SETTINGS.init_days + SETTINGS.test_days])
            for var in model.getVars():
                var_name = re.split(r'\W+', var.varName)[0]
                # if var_name == "o":
                #     index0 = int(re.split(r'\W+', var.varName)[1])
                #     o[index0] = var.X
                if var_name == "a":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    a[index0, index1] = var.X
                if var_name == "b":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    b[index0, index1] = var.X

            df.loc[(day,hospitals[0].name),"nvars"] = sum([np.product(var.shape) for var in [x, y, z, a, b]])

            return df, x, y, z, a, b

        else:
            df.loc[(day,hospitals[0].name),"nvars"] = max(df.loc[(day,hospitals[0].name),"nvars"], sum([np.product(var.shape) for var in [x, y, z]]))

            return df, x, y, z

    # Write the values found to local files.
    # with open(SETTINGS.generate_filename("results") + f"x_{SETTINGS.strategy}_{hospital.htype[:3]}_{episode}-{day}.pickle", "wb") as f:
    #     pickle.dump(x, f)
    # with open(SETTINGS.generate_filename("results") + f"y_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(y, f)
    # with open(SETTINGS.generate_filename("results") + f"z_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(z, f)

    # stop = time.perf_counter()
    # print(f"loading model results: {(stop - start):0.4f} seconds")


def model_output_to_transports(model, size_I, size_H):

    # Get the values of the model variable as found for the optimal solution.
    x = np.zeros([size_I, size_H])
    for var in model.getVars():
        var_name = re.split(r'\W+', var.varName)[0]
        if var_name == "x":
            index0 = int(re.split(r'\W+', var.varName)[1])
            index1 = int(re.split(r'\W+', var.varName)[2])
            x[index0, index1] = var.X

    return x

