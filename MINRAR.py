# from inventory import *
# from allocation import *
# from lp import *
# import csv

# avg_daily_demand, demand_name, supply_name, demand_range, supply_range, supply_size, inv_size, duration, init_length
def minrar(SETTINGS):
    
    # inventory_size = inv_size
    # demand_scenario_names = [demand_name+str(i) for i in range(demand_range[0], demand_range[1]+1)]

    SETTINGS.demand_scenario_names
    SETTINGS.supply_scenario_names

    run_RL(SETTINGS)

    # print("Supply size: ", supply_size)
    # print("Inventory size: ", inventory_size)

    # supply_scenarios = [supply_scenario(supply_size, supply_name+str(i)) for i in range(supply_range[0], supply_range[1]+1)] 
    # demand_scenarios = [multi_day_demand_scenario(396, name, avg_daily_demand) for name in demand_scenario_names]

    # A3 = [ALLOIMMUNIZATION[i]>0 or i<3 for i in range(len(ANTIGENS))] #All nonzero allo antigens

    # antigens_considered_ronald = A3
    
    
    # f = open("results/multiple(python)-ILP-"+str(inventory_size)+"-"+str(avg_daily_demand)+"#"+demand_name + "#"+supply_name+"@"+str(duration)+"-"+str(init_length)+"mip=(0.0001).csv", 'w', newline='')
    # with f:
    #     writer = csv.writer(f)

    #     #header for csv
    #     header = ["Supply Scenario", "Demand Scenario", "Daily Demand", "Inventory Size", "Duration", "Init Length", "Avg Alloimmunization Ronald", "Avg Outdates Ronald", "Avg Shorts Ronald", "Avg Num Units Short Ronald","Avg Demand Ronald"] \
    #         + ["Inventory Prevalence Ronald "+a for a in ANTIGENS] \
    #         + ["ABOD Prevalence Ronald "+a for a in ABOD_NAMES] \
    #         + ["ABOD Issuing Age Ronald "+a for a in ABOD_NAMES] \
    #         + ["ABOD Outdate % Ronald "+a for a in ABOD_NAMES] \
    #         + ["ABOD Short % Ronald "+a for a in ABOD_NAMES] \
    #         + ["Matrix Ronald "+i+"->"+j for i in ABOD_NAMES for j in ABOD_NAMES] \
    #         + ["Hist Ronald "+ str(i/400) for i in range(400)] \
    #         + ["Avg Allo per num units Ronald" + str(i+1) for i in range(4)] \
    #         + ["Num " + e + " Patients Ronald" for e in ETHNICITIES] \
    #         + ["Avg Allo " + e + " Ronald" for e in ETHNICITIES] \
    #         + ["From " + e_from + " to " + e_to + " Ronald" for e_from in ETHNICITIES for e_to in ETHNICITIES] \
    #         + ["Mismatch " + e + " " + a + " Ronald" for e in ETHNICITIES for a in ANTIGENS]

    #     writer.writerow(header)

    #     for supply_s in supply_scenarios:
    #         #load supply scenario into memory
    #         supply_s.read_from_csv()

    #         for demand_scenario in demand_scenarios:
    #             #load demand scenario into memory
    #             demand_scenario.read_from_csv()

    #             #check if there is enough supply for demand
    #             total_demand = sum(demand_scenario.demand[demand_scenario.demand["Day Needed"]<duration]["Num Units"])
    #             if (total_demand > supply_s.size):
    #                 print("ERROR, skipping "+demand_scenario.name+" because its demand "+str(total_demand)+" surpasses the supply limit of "+str(supply_s.size))
    #                 continue
                
    #             #simulation ronald
    #             print("Start Ronald", supply_s.name, demand_scenario.name)

    #             s = simulation(inventory_size, supply_s, demand_scenario, "Ronald")
    #             allo_r, outdates_r, shortage_r, num_units_short_r, avg_daily_demand_r, normalized_antigens_r, normalized_abod_r, normalized_abod_age_r, abod_outdate_r, abod_short_r, matrix_r, hist_r,allo_per_num_r,num_ethnicity,allo_per_ethnicity,ethnicity_mat, ethnicity_mismatch_mat = s.run(duration, init_length, antigens_considered_ronald)

    #             #write results to csv
    #             row = [supply_s.name, demand_scenario.name, avg_daily_demand, inventory_size, duration, init_length, allo_r, outdates_r, shortage_r, num_units_short_r, avg_daily_demand_r]+normalized_antigens_r+normalized_abod_r+normalized_abod_age_r+abod_outdate_r+abod_short_r+matrix_r+hist_r+allo_per_num_r+num_ethnicity+allo_per_ethnicity+ethnicity_mat+ethnicity_mismatch_mat
    #             writer.writerow(row)


"""Simulation for simulating online requests for products"""
class simulation():
    def __init__(self, inventory_size, supply_scenario, demand_scenario, variant):
        self.inventory_size = inventory_size
        self.supply_scenario = supply_scenario
        self.demand_scenario = demand_scenario
        self.inventory = inventory(self.inventory_size)
        self.inventory.set_scenarios(supply_scenario, demand_scenario)
        self.variant = variant


    def run(self, duration, init_length, antigens_considered, start_age=0):
        self.antigens_considered = antigens_considered

        total_demand = 0

        #counter of days
        daycount = 0

        #arrays for keeping results
        outdates = []
        shortages = []
        alloimmunization = []
        number_of_products_matched_abod = [0]*8
        number_of_patients_satisfied = 0
        total_number_of_patients = 0
        total_number_of_products_supplied = self.inventory_size
        total_units_short = 0
        antigen_presence = [0]*len(ANTIGENS)
        abod_precence = [0]*8
        abod_ages = [0]*8
        num_abod_products = [0]*8
        num_abod_requested = [0]*8
        self.abod_issuing_matrix = [[0]*8 for i in range(8)]
        self.allo_hist = [0]*400
        self.avg_allo_per_num_units_demanded = [0]*4
        self.num_patients_per_num_units_demanded = [0]*4

        self.num_patients_per_ethnicity = [0,0,0]
        self.allo_per_ethnicity = [0,0,0]
        self.ethnicity_issuing_matrix = [[0]*3 for i in range(3)]
        self.ethnicity_mismatch_matrix = [[0]*17 for i in range(3)]


        #initialize inventory with products, all with age start_age
        self.inventory.init(start_age)

        #run simulation
        while(daycount < duration):
            print(daycount, "/", duration)
            #read requets
            requests = self.inventory.read_requests_from_demand_scenario(daycount)

            #compute optimal allocation
            out, short, num_units_short, allo, units_issued = self.allocate(requests)
            #log results
            if daycount>=init_length:
                total_demand += sum(r.amount for r in requests)

                outdates.append(out)
                shortages.append(short)
                alloimmunization.append(allo)

                number_of_patients_satisfied += len(requests)-len(short)
                total_number_of_patients += len(requests)
                total_number_of_products_supplied += len(out)+len(units_issued) #add number of replenishments to total_supplied
                total_units_short += num_units_short

                for unit in self.inventory.units:
                    #log antigen precense
                    for a in range(len(ANTIGENS)):
                        antigen_presence[a] += unit.blood.vector[a]

                    #log abod distribution
                    abod_precence[unit.blood.major_index] += 1

                for unit in units_issued:
                    #log age of issuing
                    abod_ages[unit.blood.major_index] += unit.age + 1
                    number_of_products_matched_abod[unit.blood.major_index] += 1

                    #log number of products going through inventory by ABOD (outdated or issued)
                    num_abod_products[unit.blood.major_index] += 1
                
                for unit in out:
                    #log number of products going through inventory by ABOD (outdated or issued)
                    num_abod_products[unit.blood.major_index] += 1

                for request in requests:
                    #log number of units requested, categorized by ABOD
                    num_abod_requested[request.blood.major_index] += 1

            daycount += 1


        #compute aggregates
        total_allo = sum(alloimmunization)

        total_short = sum(len(i) for i in shortages)
        total_outdate = sum(len(i) for i in outdates)

        effective_duration = duration-init_length

        #compute averages
        avg_allo = total_allo/number_of_patients_satisfied
        avg_short = total_short/total_number_of_patients
        avg_outdate = total_outdate/total_number_of_products_supplied
        avg_num_units_short = total_units_short/effective_duration


        total_num_units_issued = sum(number_of_products_matched_abod)


        normalized_antigens = [round(i/(self.inventory_size*effective_duration),3) for i in antigen_presence]

        normalized_abod = [round(i/(self.inventory_size*effective_duration),3) for i in abod_precence]
        normalized_abod_age = [round(abod_ages[i]/max(number_of_products_matched_abod[i],1),3) for i in range(8)]

         #compute num oudates per major
        abod_outdate = [0]*8
        for out in outdates:
            for unit in out:
                abod_outdate[unit.blood.major_index] += 1
        num_abod_products = [max(1,i) for i in num_abod_products]
        abod_outdate = [abod_outdate[i] / num_abod_products[i] for i in range(8)]

        #compute num shortages per requested major
        abod_short = [0]*8
        for short in shortages:
            for blood in short:
                abod_short[blood.blood.major_index] += 1
        num_abod_requested = [max(1,i) for i in num_abod_requested]
        abod_short = [abod_short[i] / num_abod_requested[i] for i in range(8)]

        #normalize abod_issuing matrix
        for i in range(8):
            for j in range(8):
                self.abod_issuing_matrix[i][j] /= total_num_units_issued 

        #flatten matrix
        matrix = self.abod_issuing_matrix[0] + self.abod_issuing_matrix[1] + self.abod_issuing_matrix[2] + self.abod_issuing_matrix[3] + self.abod_issuing_matrix[4] \
            + self.abod_issuing_matrix[5] + self.abod_issuing_matrix[6] + self.abod_issuing_matrix[7]

        total_satisfied_patients = sum(self.allo_hist)

        #normalize hist
        for i in range(400):
            self.allo_hist[i] /= total_satisfied_patients



        #normalize avg_allo_per_num_units_demanded
        for i in range(4):
            self.avg_allo_per_num_units_demanded[i] /= max(self.num_patients_per_num_units_demanded[i],1)

        #compute average demand over effective duration
        d = self.demand_scenario.demand[self.demand_scenario.demand["Day Needed"] < duration]
        d = d[d["Day Needed"] >= init_length]
        avg_daily_demand = sum(d["Num Units"])/effective_duration
        
        self.allo_per_ethnicity = [self.allo_per_ethnicity[i]/max(self.num_patients_per_ethnicity[i],1) for i in range(len(self.allo_per_ethnicity))]
        
        ethnicity_mat = self.ethnicity_issuing_matrix[0]+self.ethnicity_issuing_matrix[1]+self.ethnicity_issuing_matrix[2]
        ethnicity_mismatch_mat = self.ethnicity_mismatch_matrix[0]+self.ethnicity_mismatch_matrix[1]+self.ethnicity_mismatch_matrix[2]
        return (avg_allo, avg_outdate, avg_short, avg_num_units_short, avg_daily_demand, normalized_antigens, normalized_abod,normalized_abod_age, abod_outdate, abod_short, matrix, self.allo_hist, self.avg_allo_per_num_units_demanded, self.num_patients_per_ethnicity, self.allo_per_ethnicity, ethnicity_mat,ethnicity_mismatch_mat)

    """Solve daily allocation problem and manage inventory after"""
    def allocate(self,requests):
        shortage_cost = (len(requests))+1

        if (self.variant == "Ronald"):
            assignments, shorts = solve_online_fast(requests, self.inventory.units, fifo_exponential_abod_usab, self.antigens_considered, short_weight=shortage_cost, mismatch_terms="direct-and-unnecessary", ethnic_dist=self.demand_scenario.patient_ethnic_distribution)
            num_units_short = 0
        elif (self.variant == "Joost"):
            assignments, shorts, num_units_short = solve_online_joost(requests, self.inventory.units, self.antigens_considered, short_weight=shortage_cost)

        units_issued = [self.inventory.units[j] for i,j in assignments]

        allo = self.compute_alloimmunization(assignments, requests, self.inventory.units)

        outdate,shortage = self.update_inventory(assignments, requests)

        return (outdate,shortage, num_units_short,allo,units_issued)

    """compute real alloimmunization, also for antigens ignored"""
    def compute_alloimmunization(self, assignments, requests, units):
        total_allo = 0
        for i in range(len(requests)):
            blood_vec = requests[i].blood.vector

            #units assinged to request i
            units_assigned = [j for _i,j in assignments if i==_i] 

            #if any mismatch, count allo once
            a = sum(int(any(units[j].blood.vector[k+3] for j in units_assigned) > blood_vec[k+3])*ALLOIMMUNIZATION[k+3] for k in range(len(ALLOIMMUNIZATION)-3)) 
            


            total_allo += a

            #hist logging
            index = min(int(a*400), 399)
            self.allo_hist[index] += 1

            number = len(units_assigned)-1

            request_ethnicity_index = requests[i].blood.ethnicity_index

            if (number >= 0):
                self.avg_allo_per_num_units_demanded[number] += a
                self.num_patients_per_num_units_demanded[number] += 1
                self.num_patients_per_ethnicity[request_ethnicity_index] += 1
                self.allo_per_ethnicity[request_ethnicity_index] += a
                for j in units_assigned:
                    self.ethnicity_issuing_matrix[units[j].blood.ethnicity_index][request_ethnicity_index] += 1
                for k in range(17):
                    self.ethnicity_mismatch_matrix[request_ethnicity_index][k] += int(any(units[j].blood.vector[k]>blood_vec[k] for j in units_assigned))

        return total_allo
    
    """Determine outdates and shortages and replenish inventory"""
    def update_inventory(self, pairs, requests):
        units_assigned = [j for i,j in pairs]
        units_unassigned = [j for j in range(len(self.inventory.units)) if j not in units_assigned]

        requests_assigned = [i for i,j in pairs]
        shortages = [requests[i] for i in range(len(requests)) if i not in requests_assigned]

        #log abod issuing matrix
        for i,j in pairs:
            self.abod_issuing_matrix[self.inventory.units[j].blood.major_index][requests[i].blood.major_index] += 1


        #check all unassigned units for outdating
        outdates = []
        left_over = []
        for j in units_unassigned:
            self.inventory.units[j].age += 1
            if self.inventory.units[j].age > 34:
                outdates.append(self.inventory.units[j])
            else:
                left_over.append(self.inventory.units[j])
        
        #fill inventory back up with new units
        new_required = self.inventory_size - len(left_over)
        self.inventory.units = left_over
        new_units = self.inventory.read_units_from_supply_scenario(new_required)
        self.inventory.units += new_units

        return (outdates,shortages)



