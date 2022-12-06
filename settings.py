import pandas as pd

class Settings():

    def __init__(self):

        self.home_dir = "C:/Users/Merel/Documents/Sanquin/Projects/RBC matching/Reinforcement Learning/blood_matching_RL/"

        # "demand": generate demand data
        # "supply": generate supply data
        # "optimize": run simulations and optimize matching
        self.mode = "optimize"

        #########################
        # OPTIMIZATION SETTINGS #
        #########################

        # "ILP": Use linear optimization.
        # "RL": Use reinforcement learning.
        self.method = "ILP"

         # "train" for training the RL model
        # "test" for running simulations with saved model
        self.RL_mode = "train"

        # Name of the model to be used for saving files (start with _ or leave empty).
        self.model_name = ""

        #######################
        # GENERATION SETTINGS #
        #######################

        # "regional": Use the patient group distribution of the OLVG, a regional hospital.
        # "university": Use the patient group distribution of the AMC, a university hospital.
        # "Other": Only sample patients that do not belong to a particular patient group.
        self.demand_scenario = "university"


        #########################
        # SIMULATION PARAMETERS #
        #########################

        self.test_days = 365
        self.init_days = 0

        # used for generating supply
        self.donor_eth_distr = [1, 0, 0]  # [Caucasian, African, Asian]
        self.supply_size = 40000

        # used for generating demand
        self.avg_daily_demand = 100
        self.inventory_size = 3 * self.avg_daily_demand

        # self.episodes = 1500
        self.episodes = 5

        # "major": Only match on the major antigens.
        # "relimm": Use relative immunogenicity weights for mismatching.
        # "patgroups": Use patient group specific mismatching weights.
        self.strategy = "patgroups"
        self.patgroup_musts = False

        self.demand_scenario_names = [f"regional_{i}" for i in range(5)]
        self.supply_scenario_names = [f"cau{round(self.donor_eth_distr[0]*100)}_afr{round(self.donor_eth_distr[1]*100)}_asi{round(self.donor_eth_distr[2]*100)}" for i in range(5)]


        ##########################
        # REINFORCEMENT LEARNING #
        ##########################

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
        self.gurobi_timeout = 60


    # Generate a file name for exporting log or result files.
    def generate_filename(self, output_type, model_type):

        return self.home_dir + f"{output_type}/{model_type}_"


    def create_output_file(self, PARAMS, episode):

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
        header = ["episode", "day", "name", "supply scenario", "demand scenario", "avg daily demand", "inventory size", "test days", "init days"]

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

        # Which products were issued to which patiens.
        header += ["avg issuing age"]
        header += [f"{i} to {j}" for i in ABOD_names for j in ABOD_names]
        header += [f"{eth0} to {eth1}" for eth0 in ethnicities for eth1 in ethnicities]

        # Matching performance.
        header += ["num outdates"] + [f"num outdates {i}" for i in ABOD_names]
        header += ["num shortages"] + [f"num shortages {i}" for i in ABOD_names]
        header += [f"num shortages {p}" for p in patgroups] + [f"num {p} {i+1} units short" for p in patgroups for i in range(4)]
        header += [f"num mismatches {p} {k}" for p in patgroups for k in antigens] + [f"num mismatched units {p} {k}" for p in patgroups for k in antigens]
        header += [f"num mismatches {eth} {k}" for eth in ethnicities for k in antigens]

        df = pd.DataFrame(columns = header)
        df["day"] = days
        df.index = df["day"]
        df = df.fillna(0).drop(columns=["day"])

        ##################
        # ADD BASIC INFO #
        ##################

        df.loc[days,"episode"] = episode
        df.loc[days,"name"] = self.model_name
        df.loc[days,"supply scenario"] = f"cau{round(self.donor_eth_distr[0]*100)}_afr{round(self.donor_eth_distr[1]*100)}_asi{round(self.donor_eth_distr[2]*100)}_{episode}"
        df.loc[days,"demand scenario"] = f"{self.demand_scenario}_{episode}"
        df.loc[days,"avg daily demand"] = self.avg_daily_demand
        df.loc[days,"inventory size"] = self.inventory_size
        df.loc[days,"test days"] = self.test_days
        df.loc[days,"init days"] = self.init_days

        return df


        # i = 0
        # file = f"{generate_filename("results", "ilp")}_{i}.csv"
        # while os.path.exists(file):
        #     i += 1
        #     file = f"{file[:-4]}_{i}.csv"
            
        # df.to_csv(file, sep=',', index=False)

        # self.output_index = i


