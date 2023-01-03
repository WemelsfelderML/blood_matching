import pandas as pd

class Settings():

    def __init__(self):

        self.home_dir = "C:/Users/Merel/Documents/Sanquin/Projects/RBC matching/Reinforcement Learning/blood_matching_RL/"

        # "demand": generate demand data
        # "supply": generate supply data
        # "optimize": run simulations and optimize matching
        self.mode = "optimize"
        # TODO: automatically create data before optimizing if not yet available in folders

        #########################
        # OPTIMIZATION SETTINGS #
        #########################

        # "ILP": Use linear optimization.
        # "RL": Use reinforcement learning.
        self.method = "ILP"

        # "on": online optimization.
        # "off": offline optimization.
        self.line = "on"

        #########################
        # SIMULATION PARAMETERS #
        #########################

        # self.test_days = 5 * 35
        self.test_days = 365
        self.init_days = 0

        self.episodes = (0,5)
        # (0,5): O-
        # (5,10): O+
        # (10,15): A-
        # (15,20): A+
        # (20,25): B-
        # (25,30): B+
        # (30,35): AB-
        # (35,40): AB+

        # Number of hospitals considered. If more than 1 (regional and university combined), a distribution center is included.
        # "regional": Use the patient group distribution of the OLVG, a regional hospital, with average daily demand of 50 products.
        # "university": Use the patient group distribution of the AMC, a university hospital, with average daily demand of 100 products.
        self.n_hospitals = {
            "regional" : 1,
            "university" : 0,
            "manual" : 0
        }

        self.avg_daily_demand = {
            "regional" : 50,
            "university" : 100,
            "manual" : 10
        }

        if sum(self.n_hospitals.values()) > 1:
            self.inv_size_factor = 5
        else:
            self.inv_size_factor = 3


        # Name of the model to be used for saving files.
        # self.model_name = f"I-R verhouding {self.inv_size_factor}x"
        self.model_name = "changing initial inventories"

        # "major": Only match on the major antigens.
        # "relimm": Use relative immunogenicity weights for mismatching.
        # "patgroups": Use patient group specific mismatching weights.
        self.strategy = "patgroups"
        self.patgroup_musts = True
 

        ##############################
        # GENERATING DEMAND / SUPPLY #
        ##############################

        self.donor_eth_distr = [1, 0, 0]  # [Caucasian, African, Asian]
        self.supply_size = (self.init_days + self.test_days + (self.inv_size_factor * 2)) * sum([self.n_hospitals[htype] * self.avg_daily_demand[htype] for htype in self.n_hospitals.keys()])


        ##########################
        # REINFORCEMENT LEARNING #
        ##########################

        # "train" for training the RL model
        # "test" for running simulations with saved model
        self.RL_mode = "train"

        self.nn_update_iter = 5
        self.max_memory_size = 1000
        self.train_iters = 10
        self.batch_size = 20
        self.gamma = 0.99
        self.alpha = 0.01
        self.dropout = 0.2
        self.exploration_max = 1.0
        self.exploration_min = 0.01
        self.exploration_decay = 0.1


        ####################
        # GUROBI OPTIMIZER #
        ####################

        self.show_gurobi_output = False
        self.gurobi_threads = 12
        self.gurobi_timeout = 12 * 60 * 60


    # Generate a file name for exporting log or result files.
    def generate_filename(self, output_type):

        return self.home_dir + f"{output_type}/{self.model_name}/{self.method.lower()}_"


    def initialize_output_dataframe(self, PARAMS, hospitals, episode):

        ##########
        # PARAMS #
        ##########

        antigens = PARAMS.major + PARAMS.minor
        ABOD_names = PARAMS.ABOD
        patgroups = PARAMS.patgroups
        ethnicities = ["Caucasian", "African", "Asian"]
        days = [i for i in range(self.init_days + self.test_days)]

        ##########
        # HEADER #
        ##########

        # General information.
        header = ["day", "location", "model name", "supply scenario", "demand scenario", "avg daily demand", "inventory size", "test days", "init days"]

        # Gurobi optimizer info.
        header += ["gurobi status", "nvars", "calc time"]
        header += ["objval shortages", "objval mismatches", "objval substitution", "objval fifo", "objval usability"]
        
        # Information about patients, donors, demand and supply.
        header += ["num patients"] + [f"num {eth} patients" for eth in ethnicities]
        header += [f"num {p} patients" for p in patgroups]
        header += ["num units requested"]  + [f"num units requested {p}" for p in patgroups]
        header += [f"num requests {i+1} units" for i in range(4)]
        header += ["num supplied products"] + [f"num supplied {i}" for i in ABOD_names] + [f"num requests {i}" for i in ABOD_names]
        header += [f"num {i} in inventory" for i in ABOD_names]

        if self.line == "off":
            header += ["products available today", "products in inventory today"]

        # Which products were issued to which patiens.
        header += ["avg issuing age"]
        header += [f"{i} to {j}" for i in ABOD_names for j in ABOD_names]
        header += [f"{eth0} to {eth1}" for eth0 in ethnicities for eth1 in ethnicities]
        header += [f"num allocated at dc {p}" for p in patgroups]

        # Matching performance.
        header += ["num outdates"] + [f"num outdates {i}" for i in ABOD_names]
        header += ["num shortages"] + [f"num shortages {i}" for i in ABOD_names]
        header += [f"num shortages {p}" for p in patgroups] + [f"num {p} {i+1} units short" for p in patgroups for i in range(4)] + ["num unavoidable shortages"]
        header += [f"num mismatches {p} {k}" for p in patgroups for k in antigens] + [f"num mismatched units {p} {k}" for p in patgroups for k in antigens]
        header += [f"num mismatches {eth} {k}" for eth in ethnicities for k in antigens]

        df = pd.DataFrame(columns = header)

        locations = []
        for hospital in hospitals:
            locations += [hospital.name] * len(days)

        if len(hospitals) == 1:
            df["day"] = days
        else:
            df["day"] = days * (1 + len(hospitals))
            locations += [f"dc_{episode}"] * len(days)
            
        df["location"] = locations

        df = df.set_index(['day', 'location'])
        df = df.fillna(0)

        ##################
        # ADD BASIC INFO #
        ##################

        df.loc[:,"model name"] = self.model_name
        df.loc[:,"test days"] = self.test_days
        df.loc[:,"init days"] = self.init_days
        df.loc[:,"supply scenario"] = f"cau{round(self.donor_eth_distr[0]*100)}_afr{round(self.donor_eth_distr[1]*100)}_asi{round(self.donor_eth_distr[2]*100)}_{episode}"
        
        for hospital in hospitals:
            indices = [(day,hospital.name) for day in days]
            df.loc[indices,"demand scenario"] = f"{hospital.htype}_{episode}"
            df.loc[indices,"avg daily demand"] = hospital.avg_daily_demand
            df.loc[indices,"inventory size"] = hospital.inventory_size
        
        return df


        # i = 0
        # file = f"{generate_filename("results")}_{i}.csv"
        # while os.path.exists(file):
        #     i += 1
        #     file = f"{file[:-4]}_{i}.csv"
            
        # df.to_csv(file, sep=',', index=False)

        # self.output_index = i


