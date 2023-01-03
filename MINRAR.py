from gurobipy import *
import numpy as np
import time
import math

from blood import *

class MINRAR():

    def __init__(self, SETTINGS, PARAMS):

        # All antigens mapped to their index in the list of considered antigens.
        antigens = PARAMS.major + PARAMS.minor
        self.antigens = antigens
        self.A = {antigens[k] : k for k in range(len(antigens))}   
        self.A_minor = {k : self.A[k] for k in PARAMS.minor}
        self.A_no_Fyb = {k : self.A[k] for k in antigens if k != "Fyb"}

        # All patient groups mapped to their index in the list of patient groups.
        self.P = {PARAMS.patgroups[i] : i for i in range(len(PARAMS.patgroups))}

        # Mismatching weights for each antigen.
        if "patgroups" in SETTINGS.strategy:
            self.w = np.array(PARAMS.patgroup_weights.loc[PARAMS.patgroups, antigens])
        elif "relimm" in SETTINGS.strategy:
            self.w = np.array(PARAMS.relimm_weights[antigens])[0]
        elif "major" in SETTINGS.strategy:
            self.w = np.array([0] * len(antigens))
        

    def log_results(self, SETTINGS, PARAMS, df, model, hospital, day, x=[], y=[], z=[], a=[], b=[]):

        start = time.perf_counter()

        name = hospital.name

        I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}
        R = {r : hospital.requests[r] for r in range(len(hospital.requests))}

        ABOD_names = PARAMS.ABOD
        patgroups = PARAMS.patgroups
        ethnicities = ["Caucasian", "African", "Asian"]
        r_today = [r for r in R.keys() if R[r].issuing_day == day]

        df.loc[(day,name),"num patients"] = len(r_today)
        df.loc[(day,name),"num units requested"] = sum([R[r].num_units for r in r_today])
        for eth in ethnicities:
            df.loc[(day,name),f"num {eth} patients"] = sum([1 for r in r_today if R[r].ethnicity == eth])
        for p in patgroups:
            df.loc[(day,name),f"num {p} patients"] = sum([1 for r in r_today if R[r].patgroup == p])
            df.loc[(day,name),f"num units requested {p}"] = sum([R[r].num_units for r in r_today if R[r].patgroup == p])
            df.loc[(day,name),f"num allocated at dc {p}"] = sum([R[r].num_units * R[r].allocated_from_dc for r in [r for r in R.keys() if R[r].issuing_day == (day + 1)] if R[r].patgroup == p])
        
        for u in range(1,5):
            df.loc[(day,name),f"num requests {u} units"] = sum([1 for r in r_today if R[r].num_units == u])

        df.loc[(day,name),"num supplied products"] = sum([1 for ip in I.values() if ip.age == 0])
        for major in ABOD_names:
            df.loc[(day,name),f"num supplied {major}"] = sum([1 for ip in I.values() if ip.major == major and ip.age == 0])
            df.loc[(day,name),f"num requests {major}"] = sum([1 for r in r_today if R[r].major == major])
            df.loc[(day,name),f"num {major} in inventory"] = sum([1 for ip in I.values() if ip.major == major])

        df.loc[(day,name),"objval shortages"] = sum(y[r] * (((len(R)-1) * (1 - min(1, R[r].issuing_day - day))) + 1) for r in R.keys())
        df.loc[(day,name),"objval fifo"] = sum(math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r] for i in I.keys() for r in R.keys())
        df.loc[(day,name),"objval usability"] = sum((I[i].get_usability(PARAMS, [hospital]) - R[r].get_usability(PARAMS, [hospital])) * x[i,r] for i in I.keys() for r in R.keys())
        if "patgroups" in SETTINGS.strategy:
            df.loc[(day,name),"objval mismatches"] = sum(z[r,k] * self.w[self.P[R[r].patgroup],k] for r in R.keys() for k in self.A.values())
            df.loc[(day,name),"objval substitution"] = sum(x[i,r] * self.w[self.P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k] for i in I.keys() for r in R.keys() for k in self.A_minor.values())
        else:
            df.loc[(day,name),"objval mismatches"] = sum(z[r,k] * self.w[k] for r in R.keys() for k in self.A.values())
            df.loc[(day,name),"objval substitution"] = sum(x[i,r] * self.w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for i in I.keys() for r in R.keys() for k in self.A_minor.values())

        xi = x.sum(axis=1)  # For each inventory product i∈I, xi[i] = 1 if the product is issued, 0 otherwise.
        xr = x.sum(axis=0)  # For each request r∈R, xr[r] = the number of products issued to this request.

        age_sum = 0
        issued_sum = 0
        for r in r_today:
            # Get all products from inventory that were issued to request r.
            issued = np.where(x[:,r]==1)[0]
            rq = R[r]

            mismatch = {ag:0 for ag in self.antigens}
            for ip in [I[i] for i in issued]:
                age_sum += ip.age
                issued_sum += 1
                df.loc[(day,name),f"{ip.major} to {rq.major}"] += 1
                df.loc[(day,name),f"{ip.ethnicity} to {rq.ethnicity}"] += 1
            
                # Get all antigens k where product i and request r have a mismatch.
                for ag in [self.antigens[k] for k in range(len(self.antigens)) if ip.vector[k] > rq.vector[k]]:
                    # Fy(a-b-) should only be matched on Fy(a), not on Fy(b). -> Fy(b-) only mismatch when Fy(a+)
                    if (ag != "Fyb") or (rq.vector[self.antigens.index("Fya")] == 1):    
                        mismatch[ag] = 1
                        df.loc[(day,name),[f"num mismatched units {rq.patgroup} {ag}"]] += 1

            for ag in self.antigens:
                df.loc[(day,name),[f"num mismatches {rq.patgroup} {ag}"]] += mismatch[ag]
                df.loc[(day,name),[f"num mismatches {rq.ethnicity} {ag}"]] += mismatch[ag]

        df.loc[(day,name),f"avg issuing age"] = age_sum / max(1, issued_sum)

        for ip in [I[i] for i in I.keys() if (xi[i] == 0) and (I[i].age >= (PARAMS.max_age-1))]:
            df.loc[(day,name),"num outdates"] += 1
            df.loc[(day,name),f"num outdates {ip.major}"] += 1

        df.loc[(day,name),"num unavoidable shortages"] = max(0, sum([R[r].num_units for r in r_today]) - len(I))
        for r in [r for r in r_today if y[r] == 1]:
            rq = R[r]
            df.loc[(day,name),"num shortages"] += 1
            df.loc[(day,name),f"num shortages {rq.major}"] += 1
            df.loc[(day,name),f"num shortages {rq.patgroup}"] += 1
            df.loc[(day,name),f"num {rq.patgroup} {int(rq.num_units - xr[r])} units short"] += 1

        if SETTINGS.line == "off":
            df.loc[(day,name),"products available today"] = ",".join([str(i) for i in I.keys() if a[i,day] - b[i,day] == 1])
        
        stop = time.perf_counter()
        # print(f"writing results to dataframe: {(stop - start):0.4f} seconds")

        
        return df 


    def get_mismatch_penalty(self, SETTINGS, hospital, r, z):
        
        if "patgroups" in SETTINGS.strategy:
            return sum(z[r,k] * self.w[self.P[hospital.requests[r].patgroup],k] for k in self.A.values())
        else:
            return sum(z[r,k] * self.w[k] for k in self.A.values())


    def minrar_single_hospital(self, SETTINGS, PARAMS, hospital, day):

        start = time.perf_counter()

        ################
        ## PARAMETERS ##
        ################

        model = Model(name="model")
        if SETTINGS.show_gurobi_output == False:
            model.Params.LogToConsole = 0
        model.setParam('Threads', SETTINGS.gurobi_threads)
        model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

        I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}               # Set of all inventory products.
        R = {r : hospital.requests[r] for r in range(len(hospital.requests))}                 # Set of all requests.

        bi = [I[i].get_usability(PARAMS, [hospital]) for i in I.keys()]
        br = [R[r].get_usability(PARAMS, [hospital]) for r in R.keys()]

        # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
        C = precompute_compatibility(SETTINGS, PARAMS, I, R)
        T = timewise_possible(SETTINGS, PARAMS, I, R, day)

        # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
        t = [1 - min(1, R[r].issuing_day - day) for r in R.keys()]


        ###############
        ## VARIABLES ##
        ###############

        # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
        # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
        # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.

        x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
        y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
        z = model.addVars(len(R), len(self.A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

        model.update()


        #################
        ## CONSTRAINTS ##
        #################


        # Force y[r] to 1 if not all requested units are satisfied.
        model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())

        # For each inventory product i∈I, ensure that i can not be issued more than once.
        model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
        model.addConstrs(x[i,r] <= C[i,r] * T[i,r] for i in I.keys() for r in R.keys())

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in self.A_no_Fyb.values())

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[self.A["Fyb"]] * (1 - R[r].vector[self.A["Fyb"]]) * R[r].vector[self.A["Fya"]] for i in I.keys()) <= z[r,self.A["Fyb"]] * R[r].num_units for r in R.keys())  


        ################
        ## OBJECTIVES ##
        ################

        # Assign a higher shortage penalty to requests with today as their issuing date.
        model.setObjectiveN(expr = quicksum(y[r] * (((len(R)-1) * t[r]) + 1) for r in R.keys()), index=0, priority=1, name="shortages") 
        if "patgroups" in SETTINGS.strategy:
            model.setObjectiveN(expr = quicksum(
                                            quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) +                                                  # Mismatches on minor antigens.
                                            quicksum(
                                                quicksum(x[i,r] * self.w[self.P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k] for k in self.A_minor.values())  # Minor antigen substitution.
                                                + (math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r])                                                          # FIFO penalties.
                                                + ((bi[i] - br[r]) * x[i,r])                                                                                # Product usability on major antigens.
                                            for i in I.keys()) 
                                        for r in R.keys()), index=1, priority=0, name="other")
        else:
            model.setObjectiveN(expr = quicksum(
                                            quicksum(z[r,k] * self.w[k] for k in self.A.values()) +                                                   # Mismatches on minor antigens.
                                            quicksum(
                                                quicksum(x[i,r] * self.w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for k in self.A_minor.values())   # Minor antigen substitution.
                                                + (math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r])                                          # FIFO penalties.
                                                + ((bi[i] - br[r]) * x[i,r])                                                                # Product usability on major antigens.
                                            for i in I.keys()) 
                                        for r in R.keys()), index=1, priority=0, name="other")
        
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


    def minrar_hospital_before_dc_allocation(self, SETTINGS, PARAMS, day, hospital):

        start = time.perf_counter()

        ################
        ## PARAMETERS ##
        ################

        model = Model(name="model")
        if SETTINGS.show_gurobi_output == False:
            model.Params.LogToConsole = 0
        model.setParam('Threads', SETTINGS.gurobi_threads)
        model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

        I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}               # Set of all inventory products.
        R = {r : hospital.requests[r] for r in range(len(hospital.requests))}                 # Set of all requests.

        # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
        C = precompute_compatibility(SETTINGS, PARAMS, I, R)
        T = timewise_possible(SETTINGS, PARAMS, I, R, day)

        # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
        t = [1 - min(1, R[r].issuing_day - day) for r in R.keys()]


        ###############
        ## VARIABLES ##
        ###############

        # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
        # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
        # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.

        x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
        y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
        z = model.addVars(len(R), len(self.A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

        model.update()


        #################
        ## CONSTRAINTS ##
        #################

        # Force y[r] to 1 if not all requested units are satisfied.
        model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())

        # Make sure that each request r∈R does not receive more products than the number of units requested.
        model.addConstrs(quicksum(x[i,r] for i in I.keys()) <= R[r].num_units for r in R.keys())

        # For each inventory product i∈I, ensure that i can not be issued more than once.
        model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
        model.addConstrs(x[i,r] <= C[i,r] * T[i,r] for i in I.keys() for r in R.keys())

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in self.A_no_Fyb.values())

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[self.A["Fyb"]] * (1 - R[r].vector[self.A["Fyb"]]) * R[r].vector[self.A["Fya"]] for i in I.keys()) <= z[r,self.A["Fyb"]] * R[r].num_units for r in R.keys())  


        ################
        ## OBJECTIVES ##
        ################

        # Assign a higher shortage penalty to requests with today as their issuing date.
        # TODO instead of times 4, the penalties for today's requests might be multiplied by len(R).
        model.setObjectiveN(expr = quicksum(y[r] * (((len(R)-1) * t[r]) + 1) for r in R.keys()), index=0, priority=1, name="shortages") 
        if "patgroups" in SETTINGS.strategy:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values() for r in R.keys())                             # Mismatches on minor antigens.
                                     + quicksum(x[i,r] * self.w[self.P[R[r].patgroup],k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in R.keys() if R[r].patgroup in ["Other", "Wu45"]] for k in self.A_minor.values()),  # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        else:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[k] for k in self.A.values() for r in R.keys())                                                   # Mismatches on minor antigens.
                                     + quicksum(x[i,r] * self.w[k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in R.keys() if R[r].patgroup in ["Other", "Wu45"]] for k in self.A_minor.values()),  # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        
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


    def allocate_propagated_requests_from_dc(self, SETTINGS, PARAMS, day, inventory, requests):

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
        R = {r : requests[r][0] for r in range(len(requests))}              # Set of all requests.

        # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
        C = precompute_compatibility(SETTINGS, PARAMS, I, R)
        T = timewise_possible(SETTINGS, PARAMS, I, R, day)

        # For each request r∈R, t[r] = 1 if the issuing day is tomorrow, 0 if it lies further in the future.
        t = [1 - min(1, R[r].issuing_day - (day+1)) for r in R.keys()]


        ###############
        ## VARIABLES ##
        ###############

        # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
        # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
        # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.

        x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
        y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
        z = model.addVars(len(R), len(self.A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

        model.update()


        #################
        ## CONSTRAINTS ##
        #################

        # Force y[r] to 1 if not all requested units are satisfied.
        model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())

        # Make sure that each request r∈R does not receive more products than the number of units requested.
        model.addConstrs(quicksum(x[i,r] for i in I.keys()) <= R[r].num_units for r in R.keys())

        # For each inventory product i∈I, ensure that i can not be issued more than once.
        model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
        model.addConstrs(x[i,r] <= C[i,r] * T[i,r] for i in I.keys() for r in R.keys())

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in self.A_no_Fyb.values())

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[self.A["Fyb"]] * (1 - R[r].vector[self.A["Fyb"]]) * R[r].vector[self.A["Fya"]] for i in I.keys()) <= z[r,self.A["Fyb"]] * R[r].num_units for r in R.keys())  

        # If a product is allocated from the dc's inventory, the mismatch penalty should improve w.r.t. the match made within the hospital.
        if "patgroups" in SETTINGS.strategy:
            model.addConstrs(quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) <= R[r].best_mismatch_penalty * 0.999 for r in R.keys())
        else:
            model.addConstrs(quicksum(z[r,k] * self.w[k] for k in self.A.values()) <= R[r].best_mismatch_penalty * 0.999 for r in R.keys())

        ################
        ## OBJECTIVES ##
        ################

        # Assign a higher shortage penalty to requests with tomorrow as their issuing date.
        # TODO instead of times 4, the penalties for today's requests might be multiplied by len(R).
        model.setObjectiveN(expr = quicksum(y[r] * (((len(R)-1) * t[r]) + 1) for r in R.keys()), index=0, priority=1, name="shortages") 
        if "patgroups" in SETTINGS.strategy:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values() for r in R.keys())                 # Mismatches on minor antigens.
                                     + quicksum(x[i,r] * self.w[self.P[R[r].patgroup],k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in R.keys() if R[r].patgroup =="Wu45"] for k in self.A_minor.values()),  # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        else:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[k] for k in self.A.values() for r in R.keys())                                       # Mismatches on minor antigens.
                                     + quicksum(x[i,r] * self.w[k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in R.keys() if R[r].patgroup =="Wu45"] for k in self.A_minor.values()),  # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        
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


    def mirar_hospital_after_dc_allocation(self, SETTINGS, PARAMS, day, hospital):

        start = time.perf_counter()

        ################
        ## PARAMETERS ##
        ################

        model = Model(name="model")
        if SETTINGS.show_gurobi_output == False:
            model.Params.LogToConsole = 0
        model.setParam('Threads', SETTINGS.gurobi_threads)
        model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

        I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}               # Set of all inventory products.
        R = {r : hospital.requests[r] for r in range(len(hospital.requests))}                 # Set of all requests.
        r_remaining = [r for r in R.keys() if R[r].allocated_from_dc == 0]

        bi = [I[i].get_usability(PARAMS, [hospital]) for i in I.keys()]
        br = [R[r].get_usability(PARAMS, [hospital]) for r in R.keys()]

        # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
        C = precompute_compatibility(SETTINGS, PARAMS, I, R)
        T = timewise_possible(SETTINGS, PARAMS, I, R, day)

        # For each request r∈R, t[r] = 1 if the issuing day is today, 0 if it lies in the future.
        t = [1 - min(1, R[r].issuing_day - day) for r in R.keys()]


        ###############
        ## VARIABLES ##
        ###############

        # x: For each request r∈R and inventory product i∈I, x[i,r] = 1 if r is satisfied by i, 0 otherwise.
        # y: For each request r∈R, y[r] = 1 if request r can not be fully satisfied (shortage), 0 otherwise.
        # z: For each request r∈R and antigen k∈A, z[r,k] = 1 if request r is mismatched on antigen k, 0 otherwise.

        x = model.addVars(len(I), len(R), name='x', vtype=GRB.BINARY, lb=0, ub=1)
        y = model.addVars(len(R), name='y', vtype=GRB.BINARY, lb=0, ub=1)
        z = model.addVars(len(R), len(self.A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

        model.update()


        #################
        ## CONSTRAINTS ##
        #################


        # Force y[r] to 1 if not all requested units are satisfied.
        model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in r_remaining)

        # For each inventory product i∈I, ensure that i can not be issued more than once.
        model.addConstrs(quicksum(x[i,r] for r in r_remaining) <= 1 for i in I.keys())

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # Force x[i,r] to 0 if product i∈I is outdated before request r∈R has to be issued.
        model.addConstrs(x[i,r] <= C[i,r] * T[i,r] for i in I.keys() for r in r_remaining)

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in r_remaining for k in self.A_no_Fyb.values())

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[self.A["Fyb"]] * (1 - R[r].vector[self.A["Fyb"]]) * R[r].vector[self.A["Fya"]] for i in I.keys()) <= z[r,self.A["Fyb"]] * R[r].num_units for r in r_remaining)  

        # The mismatch penalty should not be worse than already found to be possible.
        if "patgroups" in SETTINGS.strategy:
            model.addConstrs(quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) <= R[r].best_mismatch_penalty for r in r_remaining)
        else:
            model.addConstrs(quicksum(z[r,k] * self.w[k] for k in self.A.values()) <= R[r].best_mismatch_penalty for r in r_remaining)

        ################
        ## OBJECTIVES ##
        ################

        # Assign a higher shortage penalty to requests with today as their issuing date.
        # TODO instead of times 4, the penalties for today's requests might be multiplied by len(R).
        model.setObjectiveN(expr = quicksum(y[r] * (((len(R)-1) * t[r]) + 1) for r in R.keys()), index=0, priority=1, name="shortages")
        if "patgroups" in SETTINGS.strategy:
            model.setObjectiveN(expr = quicksum(
                                            quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) +                                           # Mismatches on minor antigens.
                                                quicksum(
                                                    (math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r])                                                         # FIFO penalties.
                                                    + ((bi[i] - br[r]) * x[i,r])                                                                                    # Product usability on major antigens.
                                                for i in I.keys()) 
                                            for r in r_remaining)
                                        + quicksum(x[i,r] * self.w[self.P[R[r].patgroup],k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in r_remaining if R[r].patgroup in ["Other", "Wu45"]] for k in self.A_minor.values()),  # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        else:
            model.setObjectiveN(expr = quicksum(
                                            quicksum(z[r,k] * self.w[k] for k in self.A.values()) +                                                                  # Mismatches on minor antigens.
                                                quicksum(
                                                    (math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r])                                                          # FIFO penalties.
                                                    + ((bi[i] - br[r]) * x[i,r])                                                                                     # Product usability on major antigens.
                                                for i in I.keys()) 
                                            for r in r_remaining)
                                        + quicksum(x[i,r] * self.w[k] * 2 * (1 - I[i].vector[k]) * R[r].vector[k] 
                                        for i in I.keys() for r in [r for r in r_remaining if R[r].patgroup in ["Other", "Wu45"]] for k in self.A_minor.values()),   # Minor antigen substitution.
                                        index=1, priority=0, name="other")
        
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


    def allocate_remaining_supply_from_dc(self, SETTINGS, PARAMS, day, inventory, hospitals, supply_sizes, allocations_from_dc):

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
        H = {h : hospitals[h] for h in range(len(hospitals))}               # Set of all inventory products.
        bi = [I[i].get_usability(PARAMS, hospitals, antigens=["C", "c", "E", "e", "K", "k", "Fya", "Fyb", "Jka", "Jkb"]) for i in I.keys()]


        ###############
        ## VARIABLES ##
        ###############

        # s: For each inventory product i∈I, x[i,h] = 1 product i will be shipped to hospital h.
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

        model.setObjective(expr = quicksum((math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,h])                          # FIFO penalties.
                                            + (bi[i] * x[i,h])                                                              # Product usability on major antigens.
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


    def minrar_offline(self, SETTINGS, PARAMS, hospital, days):

        start = time.perf_counter()

        ################
        ## PARAMETERS ##
        ################

        print("Creating model.")

        model = Model(name="model")
        if SETTINGS.show_gurobi_output == False:
            model.Params.LogToConsole = 0
        model.setParam('Threads', SETTINGS.gurobi_threads)
        model.setParam('TimeLimit', SETTINGS.gurobi_timeout)

        I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}               # Set of all inventory products.
        R = {r : hospital.requests[r] for r in range(len(hospital.requests))}                 # Set of all requests.
        consecutive_Is = [(i,i+1) for i in range(len(hospital.inventory)-1)]
        consecutive_days = [(day,day+1) for day in days[:-1]]
        # print("Creating compatibility matrix.")

        # I x R matrix containing a 1 if product i∈I is compatible with request r∈R, 0 otherwise.
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
        z = model.addVars(len(R), len(self.A), name='z', vtype=GRB.BINARY, lb=0, ub=1)

        print("Creating a.")
        a = model.addVars(len(I), len(days), name='a', vtype=GRB.BINARY, lb=0, ub=1)
        print("Creating b.")
        b = model.addVars(len(I), len(days), name='b', vtype=GRB.BINARY, lb=0, ub=1)

        model.update()


        #################
        ## CONSTRAINTS ##
        #################

        new_size = 0
        old_size = 0

        cumulative_requests = 0
        for day in days:
            print("Day:", day)
            r_today = [r for r in R.keys() if R[r].issuing_day == day]

            # TODO write proof for these numbers?
            i_min = hospital.inventory_size * np.floor(day / PARAMS.max_age)
            i_max = cumulative_requests + (hospital.inventory_size * (1 + np.ceil(day / PARAMS.max_age)))

            cumulative_requests += len(r_today)

            for r in r_today:
                for i in I.keys():
                    # Only include these constraints for products that are possibly present in the inventory.
                    if (i >= i_min) and (i <= i_max):

                        # A product can only be assigned to a request if it is present in inventory at the day request r is issued.
                        model.addConstr(x[i,r] <= a[i,R[r].issuing_day] - b[i,R[r].issuing_day])
                        
                        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
                        model.addConstr(x[i,r] <= C[i,r])
                        
                        new_size += 2

                    else:
                        x[i,r].ub = 0

            for i in I.keys():
                # Only include these constraints for products that are possibly present in the inventory.
                if (i >= i_min) and (i <= i_max):

                    # b[i,day] has to become 1 as soon as product i is issued.
                    model.addConstr(quicksum(x[i,r] * R[r].issuing_day for r in R.keys()) <= b[i,day] * day)
                    if day in days[PARAMS.max_age:]:
                        # b[i,day] has to become 35 days after product i has become available at the latest
                        model.addConstr(a[i,day-PARAMS.max_age] <= b[i,day])

                    new_size += 2

                else:
                    b[i,day].ub = 0



        # The products should become available in the same order as sampled.
        model.addConstrs(a[j,day] <= a[i,day] for i,j in [(i,j) for (i,j) in consecutive_Is] for day in days)
        old_size += len(consecutive_Is) * len(days)

        # b[i,day] has to become 1 as soon as product i is issued.
        # model.addConstrs(quicksum(x[i,r] * R[r].issuing_day for r in R.keys()) <= b[i,day] * day for i in I.keys() for day in days)
        # b[i,day] has to become 35 days after product i has become available at the latest
        # model.addConstrs(a[i,day-PARAMS.max_age] <= b[i,day] for i in I.keys() for day in days[PARAMS.max_age:])
        old_size += len(I) * len(days)
        old_size += len(I) * len(days[PARAMS.max_age:])

        # A product can only be assigned to a request if it is present in inventory at the day request r is issued.
        # model.addConstrs(x[i,r] <= a[i,R[r].issuing_day] - b[i,R[r].issuing_day] for i in I.keys() for r in R.keys())
        old_size += len(I) * len(R)
        
        # At every day during the simulation, the total number of products present in inventory should sum up to the hospital's inventory size.
        model.addConstrs(quicksum(a[i,day] - b[i,day] for i in I.keys()) == hospital.inventory_size for day in days)
        # print("constraints inventory total: ", len(days))
        new_size += len(days)
        old_size += len(days)

        # For each inventory product i∈I, ensure that the product is issued within the period it is available, or it is outdated.
        # TODO: only for products i that become available more than PARAMS.max_age days before the end of the simulation? 
        model.addConstrs(quicksum(x[i,r] for r in R.keys()) <= 1 for i in I.keys())
        # print("constraints issuing max once: ", len(I))
        new_size += len(I)
        old_size += len(I)

        # Force x[i,r] to 0 if a match between product i∈I and request r∈R is incompatible on antigens that are a 'must'.
        # model.addConstrs(x[i,r] <= C[i,r] for i in I.keys() for r in R.keys())
        # print("constraints compatibility: ", len(I) * len(R))
        old_size += len(I) * len(R)

        # Force y[r] to 1 if not all requested units are satisfied.
        model.addConstrs(R[r].num_units - quicksum(x[i,r] for i in I.keys()) <= R[r].num_units * y[r] for r in R.keys())
        # print("constraints shortages: ", len(R))
        new_size += len(R)
        old_size += len(R)

        # Force z[r,k] to 1 if at least one of the products i∈I that are issued to request r∈R mismatches on antigen k∈A.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[k] * (1 - R[r].vector[k]) for i in I.keys()) <= z[r,k] * R[r].num_units for r in R.keys() for k in self.A_no_Fyb.values())
        # print("constraints mismatches no Fyb: ", len(R) * len(self.A_no_Fyb))
        new_size += len(R) * len(self.A_no_Fyb)
        old_size += len(R) * len(self.A_no_Fyb)

        # A request can only be mismatched on Fyb if it is positive for Fya.
        model.addConstrs(quicksum(x[i,r] * I[i].vector[self.A["Fyb"]] * (1 - R[r].vector[self.A["Fyb"]]) * R[r].vector[self.A["Fya"]] for i in I.keys()) <= z[r,self.A["Fyb"]] * R[r].num_units for r in R.keys())  
        # print("constraints mismatches Fyb: ", len(R))
        new_size += len(R)
        old_size += len(R)

        print()
        print("new size:", new_size)
        print("old size:", old_size)


        ################
        ## OBJECTIVES ##
        ################

        # Assign a higher shortage penalty to requests with today as their issuing date.
        model.setObjectiveN(expr = quicksum(y[r] for r in R.keys()), index=0, priority=1, name="shortages") 
        if "patgroups" in SETTINGS.strategy:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values() for r in R.keys())                  # Mismatches on minor antigens.
                                       + quicksum(1 - quicksum(x[i,r] for r in R.keys()) for i in I.keys()) * len(R)
                                       , index=1, priority=0, name="other")  # Number of outdates.   TODO: times len(R) ???
        else:
            model.setObjectiveN(expr = quicksum(z[r,k] * self.w[k] for k in self.A.values() for r in R.keys())                                        # Mismatches on minor antigens.
                                       + quicksum(1 - quicksum(x[i,r] for r in R.keys()) for i in I.keys()) * len(R)
                                       , index=1, priority=0, name="other")  # Number of outdates.   TODO: times len(R) ???
        
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

