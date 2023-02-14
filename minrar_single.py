from gurobipy import *
import numpy as np
import time
import math

from blood import *
from log import *

# Single-hospital setup: MINRAR model for matching within a single hospital.
def minrar_single_hospital(SETTINGS, PARAMS, hospital, day, df):

    start = time.perf_counter()

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

    # Get the usability of all inventory products and patient requests with respect to the
    # distribution of major blood types in the patient population.
    bi = [I[i].get_usability(PARAMS, [hospital]) for i in I.keys()]
    br = [R[r].get_usability(PARAMS, [hospital]) for r in R.keys()]

    # print([f"{I[i].index}:{bi[i]}" for i in I.keys()])
    # print(br)

    # Matrices containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
    C = precompute_compatibility(SETTINGS, PARAMS, I, R)    # The product is compatible with the request on major and manditory antigens
    T = timewise_possible(SETTINGS, PARAMS, I, R, day)      # The product is not outdated before issuing date of request.

    # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
    t = [1 - min(1, R[r].day_issuing - day) for r in R.keys()]

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

    model.Params.PoolSearchMode = 2
    model.Params.PoolSolutions = 50
    model.Params.PoolGap = 0


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
    model.ModelSense = GRB.MINIMIZE

    # CHANGE (no doubt)
    # Remove variable x[i,r] if the match is not timewise or antigen compatible.
    for r in R.keys():
        for i in I.keys():
            if (C[i,r] == 0) or (T[i,r] == 0):
                model.remove(x[i,r])
        for k in A.values():
            if R[r].vector[k] == 1:
                model.remove(z[r,k])

    #################
    ## CONSTRAINTS ##
    #################

    # Force y[r] to 1 if not all requested units are satisfied.
    # model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())
    model.addConstrs((y[r] * R[r].num_units) + quicksum(x[i,r] for i in I.keys()) >= R[r].num_units for r in R.keys())

    # For each inventory product i∈I, ensure that i can not be issued more than once.
    model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())

    # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
    # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
    # model.addConstrs(x[i,r] <= C[i,r] * T[i,r] for i in I.keys() for r in R.keys())

    # CHANGE (no doubt)
    # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
    # model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in A.values())
    model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in A_no_Fyb.values())
    model.addConstrs(quicksum(x[i,r] * I[i].vector[A["Fyb"]] * (1 - R[r].vector[A["Fyb"]]) * R[r].vector[A["Fya"]] for i in I.keys()) <= z[r,A["Fyb"]] * R[r].num_units for r in R.keys())  

    ################
    ## OBJECTIVES ##
    ################

    # CHANGE (might change)
    # Assign a higher shortage penalty to requests with today as their issuing date.
    model.setObjective(expr = quicksum(y[r] * ((len(R) * t[r]) + 1) for r in R.keys())) 
    # model.setObjective(expr = quicksum(y[r] * ((len(R) * t[r]) + 1) for r in R.keys()))
    
    # # CHANGE (might change)
    if "patgroups" in SETTINGS.strategy:
        model.setObjectiveN(expr = 5 * quicksum(z[r,k] * w[P[R[r].patgroup],k] for k in A.values() for r in R.keys())
                                    # + quicksum(math.exp(-4.852 * (PARAMS.max_age - I[i].age - 1) / PARAMS.max_age) * x[i,r] for i in I.keys() for r in R.keys())
                                    + quicksum(0.5 ** ((PARAMS.max_age - I[i].age - 1) / 5) * x[i,r] for i in I.keys() for r in R.keys())
                                    + quicksum((bi[i] - br[r]) * x[i,r] for i in I.keys() for r in R.keys())
                                    + quicksum(x[i,r] * w[P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k]
                                    for k in A_minor.values() for r in [r for r in R.keys() if R[r].patgroup in ["Wu45", "Other"]] for i in I.keys())
                                    , index=1, priority=0, name="other")
    else:
        model.setObjectiveN(expr = 5 * quicksum(z[r,k] * w[k] for k in A.values() for r in R.keys())
                                    + quicksum(0.5 ** ((PARAMS.max_age - I[i].age - 1) / 5) * x[i,r] for i in I.keys() for r in R.keys())
                                    + quicksum((bi[i] - br[r]) * x[i,r] for i in I.keys() for r in R.keys())
                                    + quicksum(x[i,r] * w[k] * (1 - I[i].vector[k]) * R[r].vector[k]
                                    for k in A_minor.values() for r in [r for r in R.keys() if R[r].patgroup in ["Wu45", "Other"]] for i in I.keys())
                                    , index=1, priority=0, name="other")

    stop = time.perf_counter()
    print(f"model initialization: {(stop - start):0.4f} seconds")

    start = time.perf_counter()
    model.optimize()
    stop = time.perf_counter()
    print(f"Optimize: {(stop - start):0.4f} seconds")

    sc = model.SolCount
    print(f"Solutions found: {sc}")

    x = np.zeros([sc, len(I), len(R)])
    y = np.zeros([sc, len(R)])
    z = np.zeros([sc, len(R), len(A)])

    for s in range(sc):

        model.Params.SolutionNumber = s
    
        for var in model.getVars():
            var_name = re.split(r'\W+', var.varName)[0]
            if var_name == "x":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                x[s, index0, index1] = var.Xn
            if var_name == "y":
                index0 = int(re.split(r'\W+', var.varName)[1])
                y[s, index0] = var.Xn
            if var_name == "z":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                z[s, index0, index1] = var.Xn

    if sc > 1:

        R_today = [r for r in R.keys() if R[r].day_issuing == day]

        # For each solution found, the mismatch penalty for requests that need to be issued today.
        if "patgroups"in SETTINGS.strategy:
            mismatch_today = {s : sum([z[s,r,k] * w[P[R[r].patgroup],k] for k in A.values() for r in R_today]) for s in range(sc)}
        else:
            mismatch_today = {s : sum([z[s,r,k] * w[k] for k in A.values() for r in R_today]) for s in range(sc)}
        best = [s for s in mismatch_today.keys() if mismatch_today[s] == min(mismatch_today.values())]
        print(best)
        
        if len(best) > 1:
            avg_age_today = {s : sum([sum([x[s,i,r] for r in R_today]) * I[i].age for i in I.keys()]) / sum(sum(x[s])) for s in best}
            best = [s for s in avg_age_today.keys() if avg_age_today[s] == max(avg_age_today.values())]
            print(best)

            if len(best) > 1:
                usab_today = {s : sum([(bi[i] - br[r]) * x[s,i,r] for i in I.keys() for r in R_today]) for s in best}
                best = [s for s in usab_today.keys() if usab_today[s] == min(usab_today.values())]
                print(best)

                if len(best) > 1:
                    if "patgroups"in SETTINGS.strategy:
                        substitution_today = {s : sum([x[s,i,r] * w[P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k] for k in A.values() for r in R_today for i in I.keys()]) for s in best}
                    else:
                        substitution_today = {s : sum([x[s,i,r] * w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for k in A.values() for r in R_today for i in I.keys()]) for s in best}
                    best = [s for s in substitution_today.keys() if substitution_today[s] == min(substitution_today.values())]
                    print(best)

        x = x[best[0]]
        y = y[best[0]]
        z = z[best[0]]

    else:

        x = x[0]
        y = y[0]
        z = z[0]

    # inventory = [I[i].index for i in I.keys()]
    # print("in inventory:",inventory) 

    # issued = []
    # for i in I.keys():
    #     for r in R.keys():
    #         if (R[r].day_issuing == day) and (x[i,r]==1):
    #             issued.append(I[i].index)
    # print("issued:",issued)

    df.loc[(day,hospital.name),"gurobi status"] = model.status
    df.loc[(day,hospital.name),"nvars"] = len(model.getVars())


    return df, x, y, z