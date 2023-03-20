import time
import numpy as np

# After obtaining the optimal variable values from the solved model, write corresponding results to a csv file.
def log_results(SETTINGS, PARAMS, df, hospital, day, x=[], y=[], z=[], a=[], b=[]):

    antigens = PARAMS.major + PARAMS.minor

    # Name of the hospital (e.g. "reg_2" or "uni_0").
    name = hospital.name

    # Gather some parameters.
    I = {i : hospital.inventory[i] for i in range(len(hospital.inventory))}
    R = {r : hospital.requests[r] for r in range(len(hospital.requests))}
    ABOD_names = PARAMS.ABOD
    patgroups = PARAMS.patgroups
    ethnicities = ["Caucasian", "African", "Asian"]

    # Most results will be calculated only considering the requests that are issued today.
    r_today = [r for r in R.keys() if R[r].day_issuing == day]

    df.loc[(day,name),"logged"] = True
    df.loc[(day,name),"num patients"] = len(r_today)                                                                            # number of patients
    df.loc[(day,name),"num units requested"] = sum([R[r].num_units for r in r_today])                                           # number of units requested
    for eth in ethnicities:
        df.loc[(day,name),f"num {eth} patients"] = sum([1 for r in r_today if R[r].ethnicity == eth])                           # number of patients per ethnicity
    for p in patgroups:
        df.loc[(day,name),f"num {p} patients"] = sum([1 for r in r_today if R[r].patgroup == p])                                # number of patients per patient group
        df.loc[(day,name),f"num units requested {p}"] = sum([R[r].num_units for r in r_today if R[r].patgroup == p])            # number of units requested per patient group
        df.loc[(day,name),f"num allocated at dc {p}"] = sum([R[r].allocated_from_dc for r in R.keys() if R[r].patgroup == p])   # number of products allocated from the distribution center per patient group

    for u in range(1,5):
        df.loc[(day,name),f"num requests {u} units"] = sum([1 for r in r_today if R[r].num_units == u])                         # number of requests asking for [1-4] units

    df.loc[(day,name),"num supplied products"] = sum([1 for ip in I.values() if ip.age == 0])                                   # number of products added to the inventory at the end of the previous day
    for major in ABOD_names:
        df.loc[(day,name),f"num supplied {major}"] = sum([1 for ip in I.values() if ip.major == major and ip.age == 0])         # number of products per major blood group added to the inventory at the end of the previous day
        df.loc[(day,name),f"num requests {major}"] = sum([1 for r in r_today if R[r].major == major])                           # number of patients per major blood group
        df.loc[(day,name),f"num {major} in inventory"] = sum([1 for ip in I.values() if ip.major == major])                     # number of products in inventory per major blood group

    # print("Objective:",sum(y[r] * ((1 - min(1, R[r].day_issuing - day)) + 1) for r in R.keys()) + sum(sum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) for r in R.keys()))
    # print("Shortages:", sum(y[r] * ((1 - min(1, R[r].day_issuing - day)) + 1) for r in R.keys()))
    # print("Mismatches:", sum(sum(z[r,k] * self.w[self.P[R[r].patgroup],k] for k in self.A.values()) for r in R.keys()))
    # print("FIFO:", sum(sum(-1 * math.exp(-4.852 * (PARAMS.max_age - I[i].age - 1) / PARAMS.max_age) * x[i,r]  for i in I.keys()) for r in R.keys()))
    # print("Usability:", sum(sum((I[i].get_usability(PARAMS, [hospital]) - R[r].get_usability(PARAMS, [hospital])) * x[i,r] for i in I.keys()) for r in R.keys()))
    # print("MAS:", sum(x[i,r] * self.w[self.P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k] for r in [r for r in R.keys() if R[r].patgroup in ["Wu45", "Other"]] for i in I.keys() for k in self.A_minor.values())))
                                

    # # Values for all individual parts of the MINRAR's objective functions.
    # df.loc[(day,name),"objval shortages"] = sum(y[r] * (((len(R)-1) * (1 - min(1, R[r].day_issuing - day))) + 1) for r in R.keys())
    # # df.loc[(day,name),"objval fifo"] = sum(sum(math.exp(-4.852 * I[i].age / PARAMS.max_age) * x[i,r] for i in I.keys()) for r in R.keys())
    # df.loc[(day,name),"objval fifo"] = sum(sum(-1 * math.exp(-4.852 * (PARAMS.max_age - I[i].age - 1) / PARAMS.max_age) * x[i,r] for i in I.keys()) for r in R.keys())
    # df.loc[(day,name),"objval usability"] = sum((I[i].get_usability(PARAMS, [hospital]) - R[r].get_usability(PARAMS, [hospital])) * x[i,r] for i in I.keys() for r in R.keys())
    # if "patgroups" in SETTINGS.strategy:
    #     df.loc[(day,name),"objval mismatches"] = sum(z[r,k] * self.w[self.P[R[r].patgroup],k] for r in R.keys() for k in self.A.values())
    #     df.loc[(day,name),"objval substitution"] = sum(x[i,r] * self.w[self.P[R[r].patgroup],k] * (1 - I[i].vector[k]) * R[r].vector[k] for i in I.keys() for r in R.keys() for k in self.A_minor.values())
    # else:
    #     df.loc[(day,name),"objval mismatches"] = sum(z[r,k] * self.w[k] for r in R.keys() for k in self.A.values())
    #     df.loc[(day,name),"objval substitution"] = sum(x[i,r] * self.w[k] * (1 - I[i].vector[k]) * R[r].vector[k] for i in I.keys() for r in R.keys() for k in self.A_minor.values())

    xi = x.sum(axis=1)  # For each inventory product i∈I, xi[i] = 1 if the product is issued, 0 otherwise.
    xr = x.sum(axis=0)  # For each request r∈R, xr[r] = the number of products issued to this request.

    # for i in I.keys():
    #     for r in R.keys():
    #         if (x[i,r] > 0) and (x[i,r] < 1):
    #             print(i, r, x[i,r])

    age_sum = 0
    issued_sum = 0
    for r in r_today:
        # Get all products from inventory that were issued to request r.
        issued = np.where(x[:,r]==1)[0]
        rq = R[r]

        mismatch = {ag:0 for ag in antigens}
        for ip in [I[i] for i in issued]:
            age_sum += ip.age
            issued_sum += 1
            df.loc[(day,name),f"{ip.major} to {rq.major}"] += 1                                 # number of products per major blood group issued to requests per major blood group
            df.loc[(day,name),f"{ip.ethnicity} to {rq.ethnicity}"] += 1                         # number of products per ethnicity issued to requests per ethnicity
            
            # Get all antigens k on which product ip and request rq are mismatched.
            for ag in [antigens[k] for k in range(len(antigens)) if ip.vector[k] > rq.vector[k]]:
                # Fy(a-b-) should only be matched on Fy(a), not on Fy(b). -> Fy(b-) only mismatch when Fy(a+)
                if (ag != "Fyb") or (rq.vector[antigens.index("Fya")] == 1):
                    mismatch[ag] = 1
                    df.loc[(day,name),[f"num mismatched units {rq.patgroup} {ag}"]] += 1        # number of mismatched units per patient group and antigen

        for ag in antigens:
            df.loc[(day,name),[f"num mismatches {rq.patgroup} {ag}"]] += mismatch[ag]           # number of mismatched patients per patient group and antigen
            df.loc[(day,name),[f"num mismatches {rq.ethnicity} {ag}"]] += mismatch[ag]          # number of mismatched patients per patient ethnicity and antigen

    df.loc[(day,name),f"avg issuing age"] = age_sum / max(1, issued_sum)                        # average age of all issued products

    for ip in [I[i] for i in I.keys() if (xi[i] == 0) and (I[i].age >= (PARAMS.max_age-1))]:
        df.loc[(day,name),"num outdates"] += 1                                                  # number of outdated inventory products
        df.loc[(day,name),f"num outdates {ip.major}"] += 1                                      # number of outdated inventory products per major blood group

    df.loc[(day,name),"num unavoidable shortages"] = max(0, sum([R[r].num_units for r in r_today]) - len(I))                # difference between the number of requested units and number of products in inventory, in case the former is larger
    for r in [r for r in r_today if y[r] == 1]:
        rq = R[r]
        df.loc[(day,name),"num shortages"] += 1                                                 # number of today's requests that were left unsatisfied
        df.loc[(day,name),f"num shortages {rq.major}"] += 1                                     # number of unsatisfied requests per major blood group
        df.loc[(day,name),f"num shortages {rq.patgroup}"] += 1                                  # number of unsatisfied requests per patient group
        df.loc[(day,name),f"num {rq.patgroup} {int(rq.num_units - xr[r])} units short"] += 1    # difference between the number units requested and issued

    if SETTINGS.line == "off":
        df.loc[(day,name),"products available today"] = ",".join([str(i) for i in I.keys() if a[i,day] - b[i,day] == 1])    # this number should be equal to the inventory size provided in the settings
    

    # Write the values found to pickle files.
    # with open(SETTINGS.generate_filename("results") + f"x_{SETTINGS.strategy}_{hospital.htype[:3]}_{episode}-{day}.pickle", "wb") as f:
    #     pickle.dump(x, f)
    # with open(SETTINGS.generate_filename("results") + f"y_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(y, f)
    # with open(SETTINGS.generate_filename("results") + f"z_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(z, f)
    
    return df 