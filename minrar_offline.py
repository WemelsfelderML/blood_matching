from gurobipy import *
import numpy as np
import time
import math

from blood import *
from log import *

# Single-hospital scenario with Offline solving.
def minrar_offline(SETTINGS, PARAMS, hospital, days):

    start = time.perf_counter()

    ############
    ## GUROBI ##
    ############

    model = Model(name="model")
    if SETTINGS.show_gurobi_output == False:
        model.Params.LogToConsole = 0
    if SETTINGS.gurobi_threads != None:
        model.setParam('Threads', SETTINGS.gurobi_threads)
    if SETTINGS.gurobi_timeout != None:
        model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

    ################
    ## PARAMETERS ##
    ################

    # The set of antigens is reshaped to several formats, for different uses within this class.
    antigens = PARAMS.major + PARAMS.minor
    A = {antigens[k] : k for k in range(len(antigens))}   
    A_minor = {k : A[k] for k in PARAMS.minor}
    A_no_Fyb = {k : A[k] for k in antigens if k != "Fyb"}

    # Mapping of the patient groups to their index in the list of patient groups.
    P = {PARAMS.patgroups[i] : i for i in range(len(PARAMS.patgroups))}

    # Retrieve the antigen (and patient group) weights.
    if "patgroups" in SETTINGS.strategy:
        w = np.array(PARAMS.patgroup_weights.loc[PARAMS.patgroups, antigens])
    elif "relimm" in SETTINGS.strategy:
        w = np.array(PARAMS.relimm_weights[antigens])[0]
    elif "major" in SETTINGS.strategy:
        w = np.array([0] * len(antigens))

    # Sets of all inventory products and patient requests.
    I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}
    R = {r : hospital.requests[r] for r in range(len(hospital.requests))}

    # List of pairs of inventory products that should be supplied in that order.
    consecutive_Is = [(i,i+1) for i in range(len(hospital.inventory)-1)]

    # Matrix containing a 1 if product i∈I is compatible with request r∈R on the major and manditory antigens, 0 otherwise.
    C = precompute_compatibility(SETTINGS, PARAMS, I, R)

    ###############
    ## VARIABLES ##
    ###############

    # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
    # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
    # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.
    # p: 1 if product i has already been sampled on some day, 0 if not (yet) sampled.
    # d: 1 if product i has already been issued, 0 if not (yet) issued.

    print("Creating x.")
    x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
    print("Creating y.")
    y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
    print("Creating z.")
    z = model.addVars(len(R), len(A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

    print("Creating a.")
    a = model.addVars(len(I), len(days), name='a', vtype=GRB.BINARY, lb=0, ub=1)
    print("Creating b.")
    b = model.addVars(len(I), len(days), name='b', vtype=GRB.BINARY, lb=0, ub=1)

    model.update()

    # Remove variable x[i,r] if the match is not timewise or antigen compatible.
    for r in R.keys():
        for i in I.keys():
            if C[i,r] == 0:
                model.remove(x[i,r])

    #################
    ## CONSTRAINTS ##
    #################

    nvars = (len(I)*len(R)) + len(R) + (len(R)*len(A)) + (len(I)*len(days)) + (len(I)*len(days))
    ncons = 0

    cumulative_requests = 0
    for day in days:
        print("Day:", day)
        r_today = [r for r in R.keys() if R[r].day_issuing == day]

        # TODO write proof for these numbers.
        # The minimum and maximum index i∈I that can currently be available in the inventory.
        i_min = hospital.inventory_size * np.floor(day / PARAMS.max_age)
        i_max = cumulative_requests + (hospital.inventory_size * (1 + np.ceil(day / PARAMS.max_age)))

        cumulative_requests += len(r_today)

        for r in r_today:

            for k in A_no_Fyb.values():
                if R[r].vector[k] == 0:
                    # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
                    model.addConstr(quicksum(x[i,r] for i in [i for i in I.keys() if I[i].vector[k] == 1]) <= z[r,k] * R[r].num_units)
                    ncons += 1
                else:
                    model.remove(z[r,k])
                    nvars -= 1
                  
            # A request can only be mismatched on Fyb if it is positive for Fya.  
            if (R[r].vector[A["Fyb"]] == 0) and (R[r].vector[A["Fya"]] == 1):
                # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen Fyb.
                model.addConstr(quicksum(x[i,r] for i in [i for i in I.keys() if I[i].vector[A["Fyb"]] == 1]) <= z[r,A["Fyb"]] * R[r].num_units)
                ncons += 1
            else:
                model.remove(z[r,A["Fyb"]])
                nvars -= 1

            for i in I.keys():
                # Only include these constraints for products that are possibly present in the inventory.
                if (i >= i_min) and (i <= i_max):

                    # A product can only be assigned to a request if it is present in inventory at the day request r is issued.
                    model.addConstr(x[i,r] <= a[i,R[r].day_issuing] - b[i,R[r].day_issuing])
                    ncons += 1

                else:
                    model.remove(x[i,r])
                    nvars -= 1

        for i in I.keys():
            # Only include these constraints for products that are possibly present in the inventory.
            if (i >= i_min) and (i <= i_max):

                # b[i,day] is forced to 1 as soon as product i is issued.
                model.addConstr(quicksum(x[i,r] * R[r].day_issuing for r in R.keys()) <= b[i,day] * day)
                if day in days[PARAMS.max_age:-(PARAMS.max_age-1)]:
                    # Max-age days after product i has become available, b[i,day] is forced to 1.
                    model.addConstr(a[i,day-PARAMS.max_age] <= b[i,day])

                ncons += 2

            else:
                model.remove(b[i,day])
                nvars -= 1

    # The products should become available in the same order as sampled.
    model.addConstrs(a[j,day] <= a[i,day] for i,j in [(i,j) for (i,j) in consecutive_Is] for day in days)
    ncons += len(consecutive_Is) * len(days)
    
    # At every day during the simulation, the total number of products present in inventory should sum up to the hospital's inventory size.
    model.addConstrs(quicksum(a[i,day] - b[i,day] for i in I.keys()) == hospital.inventory_size for day in days)
    # print("constraints inventory total: ", len(days))
    ncons += len(days)

    # For each inventory product i∈I, ensure that i can not be issued more than once.
    model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())
    ncons += len(I)

    # Force y[r] to 1 if not all requested units are satisfied.
    model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())
    # print("constraints shortages: ", len(R))
    ncons += len(R)

    # # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
    # model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in A_no_Fyb.values())
    # # print("constraints mismatches no Fyb: ", len(R) * len(A_no_Fyb))
    # ncons += len(R) * len(A_no_Fyb)

    # # A request can only be mismatched on Fyb if it is positive for Fya.
    # model.addConstrs(quicksum(x[i,r] * I[i].vector[A["Fyb"]] * (1 - R[r].vector[A["Fyb"]]) * R[r].vector[A["Fya"]] for i in I.keys()) <= z[r,A["Fyb"]] * R[r].num_units for r in R.keys())  
    # # print("constraints mismatches Fyb: ", len(R))
    # ncons += len(R)

    print()
    print("number of variables:", nvars)
    print("number of constraints:", ncons)


    ################
    ## OBJECTIVES ##
    ################

    # Assign a higher shortage penalty to requests with today as their issuing date.
    model.setObjectiveN(expr = quicksum(y[r] for r in R.keys()), index=0, priority=1, name="shortages") 
    if "patgroups" in SETTINGS.strategy:
        model.setObjectiveN(expr = quicksum(z[r,k] * w[P[R[r].patgroup],k] for k in A.values() for r in R.keys())    # Mismatches on minor antigens.
                                   + quicksum(1 - quicksum(x[i,r] for r in R.keys()) for i in I.keys()) * len(R)                    # Number of outdates.   TODO: times len(R) ???
                                   , index=1, priority=0, name="other")
    else:
        model.setObjectiveN(expr = quicksum(z[r,k] * w[k] for k in A.values() for r in R.keys())                          # Mismatches on minor antigens.
                                   + quicksum(1 - quicksum(x[i,r] for r in R.keys()) for i in I.keys()) * len(R)                    # Number of outdates.   TODO: times len(R) ???
                                   , index=1, priority=0, name="other")
    
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
    print(f"optimize: {(stop - start):0.4f} seconds")

    print(PARAMS.status_code[model.status])
    
    return model
