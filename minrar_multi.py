from gurobipy import *
import numpy as np
import time
import math

from blood import *
from log import *


# Multi-hospital setup: MINRAR model for matching simultaniously in multiple hospitals.
def minrar_multiple_hospitals(SETTINGS, PARAMS, dc, hospitals, day, df):

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

    # Sets of all hospitals, their inventory products and patient requests.
    H = {h : hospitals[h] for h in range(len(hospitals))}
    R = {h : {} for h in H.keys()}
    Ih = {h : {} for h in H.keys()}
    for h in H.keys():
        R[h] = {r : hospitals[h].requests[r] for r in range(len(hospitals[h].requests))}        # Requests for each hospital.
        Ih[h] = {i : hospitals[h].inventory[i] for i in range(len(hospitals[h].inventory))}     # Products in hospital inventories.
    Idc = {i : dc.inventory[i] for i in range(len(dc.inventory))}                               # Products in the distribution center's inventory.

    # Get the usability of all inventory products, in both the hospitals and the distribution center, and for all
    # patient requests with respect to the distribution of major blood types in the patient population.
    bih = {h : [Ih[h][i].get_usability(PARAMS, [hospitals[h]]) for i in Ih[h].keys()] for h in H.keys()}
    bidc = [Idc[i].get_usability(PARAMS, hospitals) for i in Idc.keys()]
    br = {h : [R[h][r].get_usability(PARAMS, [hospitals[h]]) for r in R[h].keys()] for h in H.keys()}

    # Matrices containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
    Ch = {h : precompute_compatibility(SETTINGS, PARAMS, Ih[h], R[h]) for h in H.keys()}    # The product in the hospital's inventory is compatible on major and manditory antigens.
    Cdc = {h : precompute_compatibility(SETTINGS, PARAMS, Idc, R[h]) for h in H.keys()}     # The product in the distribution center's inventory is compatible on major and manditory antigens.
    Th = {h : timewise_possible(SETTINGS, PARAMS, Ih[h], R[h], day) for h in H.keys()}      # The product in the hospital's inventory is not outdated before issuing date of request.
    Tdc = {h : timewise_possible(SETTINGS, PARAMS, Idc, R[h], day) for h in H.keys()}       # The product in the distribution center's inventory is not outdated before issuing date of request.

    # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
    t = [[1 - min(1, R[h][r].day_issuing - day) for r in R[h].keys()] for h in H.keys()]

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

    # For each hospital h∈H:
        # xh: For eachrequest r∈R[h] and product i∈Ih[h] (hospital's inventory), xh[h][i,r] = 1 if r is satisfied by i, 0 otherwise.
        # xdc: For each request r∈R[h] and i∈Idc (distribution center's inventory), xdc[h][i,r] = 1 if r is satisfied by i, 0 otherwise.
        # y: For each request r∈R[h], y[h][r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
        # z: For each request r∈R[h] and antigen k∈A, z[h][r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.
    xh = [model.addVars(len(Ih[h]), len(R[h]), name=f"xh{h}", vtype=GRB.BINARY, lb=0, ub=1) for h in H.keys()]
    xdc = [model.addVars(len(Idc), len(R[h]), name=f"xdc{h}", vtype=GRB.BINARY, lb=0, ub=1) for h in H.keys()]
    y = [model.addVars(len(R[h]), name=f"y{h}", vtype=GRB.BINARY, lb=0, ub=1) for h in H.keys()]
    z = [model.addVars(len(R[h]), len(A), name=f"z{h}", vtype=GRB.BINARY, lb=0, ub=1) for h in H.keys()]

    model.update()
    model.ModelSense = GRB.MINIMIZE

    for h in H.keys():
        for r in R[h].keys():

            # Remove variable xh[h][i,r] if the match is not timewise or antigen compatible.
            for i in Ih[h].keys():
                if (Ch[h][i,r] == 0) or (Th[h][i,r] == 0):
                    model.remove(xh[h][i,r])
            
            # Remove variable xdc[h][i,r] if the match is not timewise or antigen compatible.
            for i in Idc.keys():
                if (Cdc[h][i,r] == 0) or (Tdc[h][i,r] == 0):
                    model.remove(xdc[h][i,r])

            for k in A.values():
                if R[h][r].vector[k] == 1:
                    model.remove(z[h][r,k])

            # Remove variable xdc[h][i,r] if the issuing date of request r is today.
            if t[h][r] == 1:
                for i in Idc.keys():
                    model.remove(xdc[h][i,r])

    #################
    ## CONSTRAINTS ##
    #################

    ncons = 0

    for h in H.keys():

        Rh = R[h]
        Ihh = Ih[h]

        # Force y[r] to 1 if not all requested units are satisfied (either from the hospital's own inventory or from the dc's inventory).
        # model.addConstrs(Rh[r].num_units - quicksum(xh[h][i,r] for i in Ihh.keys()) - quicksum(xdc[h][i,r] for i in Idc.keys()) <= Rh[r].num_units * y[h][r] for r in Rh.keys())
        model.addConstrs((y[h][r] * Rh[r].num_units) + quicksum(xh[h][i,r] for i in Ihh.keys()) + quicksum(xdc[h][i,r] for i in Idc.keys()) >= Rh[r].num_units for r in Rh.keys())
        ncons += len(Rh)

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
        # model.addConstrs(xh[h][i,r] <= Ch[h][i,r] * Th[h][i,r] for i in Ihh.keys() for r in Rh.keys())
        # model.addConstrs(xdc[h][i,r] <= Cdc[h][i,r] * Tdc[h][i,r] for i in Idc.keys() for r in Rh.keys())
        # ncons += (len(Ihh) * len(Rh)) + (len(Idc) * len(Rh))

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(xh[h][i,r] * Ihh[i].vector[k] * (1 - Rh[r].vector[k]) for i in Ihh.keys()) <= z[h][r,k] * Rh[r].num_units for r in Rh.keys() for k in A_no_Fyb.values())
        model.addConstrs(quicksum(xdc[h][i,r] * Idc[i].vector[k] * (1 - Rh[r].vector[k]) for i in Idc.keys()) <= z[h][r,k] * Rh[r].num_units for r in Rh.keys() for k in A_no_Fyb.values())
        ncons += (len(Rh) * len(A_no_Fyb)) + (len(Rh) * len(A_no_Fyb))

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(xh[h][i,r] * Ihh[i].vector[A["Fyb"]] * (1 - Rh[r].vector[A["Fyb"]]) * Rh[r].vector[A["Fya"]] for i in Ihh.keys()) <= z[h][r,A["Fyb"]] * Rh[r].num_units for r in Rh.keys())
        model.addConstrs(quicksum(xdc[h][i,r] * Idc[i].vector[A["Fyb"]] * (1 - Rh[r].vector[A["Fyb"]]) * Rh[r].vector[A["Fya"]] for i in Idc.keys()) <= z[h][r,A["Fyb"]] * Rh[r].num_units for r in Rh.keys())
        ncons += len(Rh) + len(Rh.keys())

        # For each request, the number of products allocated by the hospital and DC together should not exceed the number of units requested.
        # model.addConstrs(quicksum(xh[h][i,r] for i in Ihh.keys()) + quicksum(xdc[h][i,r] for i in Idc.keys()) <= Rh[r].num_units for r in Rh.keys())
        # ncons += len(Rh)

        # Force xdc[i,r] to 0 if the issuing date of request r is today.
        # model.addConstrs(xdc[h][i,r] + t[h][r] <= 1 for i in Idc.keys() for r in Rh.keys())
        # ncons += len(Idc) * len(Rh)

        # For each inventory product i∈I, ensure that i can not be issued more than once.
        model.addConstrs(quicksum(xh[h][i,r] for r in Rh.keys()) <= 1 for i in Ihh.keys())
        ncons += len(Ihh)
    model.addConstrs(quicksum(quicksum(xdc[h][i,r] for r in R[h].keys()) for h in H.keys()) <= 1 for i in Idc.keys())
    ncons += len(Idc)

    print("ncons:",ncons)

    
    ################
    ## OBJECTIVES ##
    ################

    # Assign a higher shortage penalty to requests with today as their issuing date.
    model.setObjective(expr = quicksum(quicksum(y[h][r] * ((len(R[h]) * t[h][r]) + 1) for r in R[h].keys()) for h in H.keys()))         # Shortages.

    if "patgroups" in SETTINGS.strategy:
        model.setObjectiveN(expr = 5 * quicksum(quicksum(quicksum(z[h][r,k] * w[P[R[h][r].patgroup],k] for r in R[h].keys()) for h in H.keys()) for k in A.values())
                                    + quicksum(quicksum(0.5 ** ((PARAMS.max_age - Ih[h][i].age - 1) / 5) * xh[h][i,r] for i in Ih[h].keys() for r in R[h].keys()) for h in H.keys())
                                    + quicksum(quicksum(quicksum(0.5 ** ((PARAMS.max_age - Idc[i].age - 1) / 5) * xdc[h][i,r] for r in R[h].keys()) for h in H.keys()) for i in Idc.keys())
                                    + quicksum(quicksum((bih[h][i] - br[h][r]) * xh[h][i,r] for i in Ih[h].keys() for r in R[h].keys()) for h in H.keys())
                                    + quicksum(quicksum(quicksum((bidc[i] - br[h][r]) * xdc[h][i,r] for r in R[h].keys()) for h in H.keys()) for i in Idc.keys())
                                    + quicksum(quicksum(quicksum(xh[h][i,r] * w[P[R[h][r].patgroup],k] * (1 - Ih[h][i].vector[k]) * R[h][r].vector[k] for r in [r for r in R[h].keys() if R[h][r].patgroup in ["Wu45", "Other"]] for i in Ih[h].keys()) for h in H.keys()) for k in A_minor.values())
                                    + quicksum(quicksum(quicksum(xdc[h][i,r] * w[P[R[h][r].patgroup],k] * (1 - Idc[i].vector[k]) * R[h][r].vector[k] for r in [r for r in R[h].keys() if R[h][r].patgroup in ["Wu45", "Other"]]) for h in H.keys()) for i in Idc.keys() for k in A_minor.values())
                                    , index=1, priority=0, name="other")
    else:
        model.setObjectiveN(expr = 5 * quicksum(quicksum(quicksum(z[h][r,k] * w[k] for r in R[h].keys()) for h in H.keys()) for k in A.values())
                                    + quicksum(quicksum(0.5 ** ((PARAMS.max_age - Ih[h][i].age - 1) / 5) * xh[h][i,r] for i in Ih[h].keys() for r in R[h].keys()) for h in H.keys())
                                    + quicksum(quicksum(quicksum(0.5 ** ((PARAMS.max_age - Idc[i].age - 1) / 5) * xdc[h][i,r] for r in R[h].keys()) for h in H.keys()) for i in Idc.keys())
                                    + quicksum(quicksum((bih[h][i] - br[h][r]) * xh[h][i,r] for i in Ih[h].keys() for r in R[h].keys()) for h in H.keys())
                                    + quicksum(quicksum(quicksum((bidc[i] - br[h][r]) * xdc[h][i,r] for r in R[h].keys()) for h in H.keys()) for i in Idc.keys())
                                    + quicksum(quicksum(quicksum(xh[h][i,r] * w[k] * (1 - Ih[h][i].vector[k]) * R[h][r].vector[k] for r in [r for r in R[h].keys() if R[h][r].patgroup in ["Wu45", "Other"]] for i in Ih[h].keys()) for h in H.keys()) for k in A_minor.values())
                                    + quicksum(quicksum(quicksum(xdc[h][i,r] * w[k] * (1 - Idc[i].vector[k]) * R[h][r].vector[k] for r in [r for r in R[h].keys() if R[h][r].patgroup in ["Wu45", "Other"]]) for h in H.keys()) for i in Idc.keys() for k in A_minor.values())
                                    , index=1, priority=0, name="other")

    stop = time.perf_counter()
    print(f"model initialization: {(stop - start):0.4f} seconds")

    start = time.perf_counter()
    model.optimize()
    stop = time.perf_counter()
    print(f"Optimize: {(stop - start):0.4f} seconds")

    sc = model.SolCount
    print(f"Solutions found: {sc}")

    # Create numpy arrays filled with zeros.
    xh = [np.zeros([sc, len(Ih[h]), len(R[h])]) for h in H.keys()]
    xdc = [np.zeros([sc, len(Idc), len(R[h])]) for h in H.keys()]
    y = [np.zeros([sc, len(R[h])]) for h in H.keys()]
    z = [np.zeros([sc, len(R[h]), len(A)]) for h in H.keys()]

    for s in range(sc):

        model.Params.SolutionNumber = s
    
        for h in range(len(hospitals)):
            for var in model.getVars():
                var_name = re.split(r'\W+', var.varName)[0]
                if var_name == f"xh{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xh[h][s, index0, index1] = var.X
                if var_name == f"xdc{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xdc[h][s, index0, index1] = var.X
                if var_name == f"y{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    y[h][s, index0] = var.X
                if var_name == f"z{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    z[h][s, index0, index1] = var.X

    if sc > 1:

        R_today = {h : [r for r in R[h].keys() if R[h][r].day_issuing == day] for h in H.keys()}

        # For each solution found, the mismatch penalty for requests that need to be issued today.
        if "patgroups"in SETTINGS.strategy:
            mismatch_today = {s : sum([sum([z[h][s,r,k] * w[P[R[h][r].patgroup],k] for r in R_today[h]]) for h in H.keys() for k in A.values()]) for s in range(sc)}
        else:
            mismatch_today = {s : sum([sum([z[h][s,r,k] * w[k] for r in R_today[h]]) for h in H.keys() for k in A.values()]) for s in range(sc)}
        best = [s for s in mismatch_today.keys() if mismatch_today[s] == min(mismatch_today.values())]
        print(best)
        
        if len(best) > 1:
            avg_age_today = {s : sum([sum([sum([xh[h][s,i,r] for r in R_today[h]]) * Ih[h][i].age for i in Ih[h].keys()]) + sum([sum([xdc[h][s,i,r] for r in R_today[h]]) * Idc[i].age for i in Idc.keys()]) for h in H.keys()]) / sum([sum(sum(xh[h][s])) + sum(sum(xdc[h][s])) for h in H.keys()]) for s in best}
            best = [s for s in avg_age_today.keys() if avg_age_today[s] == max(avg_age_today.values())]
            print(best)

            if len(best) > 1:
                usab_today = {s : sum([sum([(bih[h][i] - br[h][r]) * xh[h][s,i,r] for i in Ih[h].keys() for r in R_today[h]]) + sum([(bidc[i] - br[h][r]) * xdc[h][s,i,r] for i in Idc.keys() for r in R_today[h]]) for h in H.keys()]) for s in best}
                best = [s for s in usab_today.keys() if usab_today[s] == min(usab_today.values())]
                print(best)

                if len(best) > 1:
                    if "patgroups"in SETTINGS.strategy:
                        substitution_today = {s : sum([sum([xh[h][s,i,r] * w[P[R[h][r].patgroup],k] * (1 - Ih[h][i].vector[k]) * R[h][r].vector[k] for k in A.values() for r in R_today[h] for i in Ih[h].keys()]) for h in H.keys()]) + sum([sum([xdc[h][s,i,r] * w[P[R[h][r].patgroup],k] * (1 - Idc[i].vector[k]) * R[h][r].vector[k] for k in A.values() for r in R_today[h] for i in Idc.keys()]) for h in H.keys()]) for s in best}
                    else:
                        substitution_today = {s : sum([sum([xh[h][s,i,r] * w[k] * (1 - Ih[h][i].vector[k]) * R[h][r].vector[k] for k in A.values() for r in R_today[h] for i in Ih[h].keys()]) for h in H.keys()]) + sum([sum([xdc[h][s,i,r] * w[k] * (1 - Idc[i].vector[k]) * R[h][r].vector[k] for k in A.values() for r in R_today[h] for i in Idc.keys()]) for h in H.keys()]) for s in best}
                    best = [s for s in substitution_today.keys() if substitution_today[s] == min(substitution_today.values())]
                    print(best)

        for h in H.keys():
            xh[h] = xh[h][best[0]]
            xdc[h] = xdc[h][best[0]]
            y[h] = y[h][best[0]]
            z[h] = z[h][best[0]]

    else:
        for h in H.keys():
            xh[h] = xh[h][0]
            xdc[h] = xdc[h][0]
            y[h] = y[h][0]
            z[h] = z[h][0]

    
    # df.loc[(day,hospital.name),"gurobi status"] = model.status
    # df.loc[(day,hospital.name),"nvars"] = len(model.getVars())


    return df, xh, xdc, y, z


# Multi-hospital setup: allocate products to each of the hospitals to restock them upto their maximum capacity.
def allocate_remaining_supply_from_dc(SETTINGS, PARAMS, day, inventory, hospitals, supply_sizes, allocations_from_dc):

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

    # Sets of all hospitals and of the distribution center's inventory products.
    H = {h : hospitals[h] for h in range(len(hospitals))}
    I = {i : inventory[i] for i in range(len(inventory))}
    
    # Get the usability of all inventory products with respect to the distribution of major blood types in the patient population.
    # bi = [I[i].get_usability(PARAMS, hospitals, antigens=["C", "c", "E", "e", "K", "k", "Fya", "Fyb", "Jka", "Jkb"]) for i in I.keys()]   # CHANGE (testing needed)
    bi = [I[i].get_usability(PARAMS, hospitals, antigens=PARAMS.minor) for i in I.keys()]


    ###############
    ## VARIABLES ##
    ###############

    # For each inventory product i∈I, x[i,h] = 1 if product i will be shipped to hospital h, 0 otherwise.
    x = model.addVars(len(I), len(hospitals), name='x', vtype=GRB.BINARY, lb=0, ub=1)
    model.update()


    #################
    ## CONSTRAINTS ##
    #################

    # Force x[i,h] to 1 if product i∈I was already allocated to hospital h∈H in the previous optimization.
    model.addConstrs(x[i,h] >= allocations_from_dc[i,h] for i in I.keys() for h in H.keys())

    # Make sure the number of supplied products is at least the necessary amount to restock each hospital completely.
    model.addConstrs(quicksum(x[i,h] for i in I.keys()) >= supply_sizes[h] for h in H.keys())

    # For each inventory product i∈I, ensure that i can not be allocated more than once.
    model.addConstrs(quicksum(x[i,h] for h in H.keys()) <= 1 for i in I.keys())


    ################
    ## OBJECTIVES ##
    ################

    model.setObjective(expr = quicksum((math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,h])      # FIFO penalties.
                                        + (bi[i] * x[i,h])                                          # Product usability on major antigens.
                                        for i in I.keys() for h in H.keys()))

    # Minimize the objective functions.
    model.ModelSense = GRB.MINIMIZE

    stop = time.perf_counter()
    # print(f"model initialization: {(stop - start):0.4f} seconds")


    ##############
    ## OPTIMIZE ##
    ##############

    start = time.perf_counter()
    model.optimize()
    stop = time.perf_counter()
    # print(f"optimize: {(stop - start):0.4f} seconds")

    print(PARAMS.status_code[model.status])

    return model


    
