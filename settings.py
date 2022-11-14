class Settings():

    def __init__(self):

        self.home_dir = "C:/Users/Merel/Documents/Sanquin/Projects/RBC matching/Reinforcement Learning/code/"

        # "demand": generate demand data
        # "supply": generate supply data
        # "optimize": run simulations and optimize matching
        self.mode = "optimize"

        self.test_days = 365
        self.init_days = 31

        # used for generating supply
        self.donor_eth_distr = [0.95, 0.05, 0]  # [Caucasian, African, Asian]
        self.supply_size = 200000

        # used for generating demand
        self.avg_daily_demand = 100
        self.max_inventory_size = 3 * self.avg_daily_demand


        # "regional": Use the patient group distribution of the OLVG, a regional hospital.
        # "university": Use the patient group distribution of the AMC, a university hospital.
        # "Other": Only sample patients that do not belong to a particular patient group.
        self.demand_scenario = "university"


        # Used for generating supply and demand
        # self.episodes = 1500
        self.episodes = 10


        if self.mode == "demand":
            self.name = self.demand_scenario
        elif self.mode == "supply":
            eth = self.donor_eth_distr
            self.name = f"cau{round(eth[0]*100)}_afr{round(eth[1]*100)}_asi{round(eth[2]*100)}"


        ##################################
        # ONLY USED IF MODE = "optimize" #

        # "ILP": Use linear optimization.
        # "RL": Use reinforcement learning.
        self.method = "RL"

        # "train" for training the RL model
        # "test" for running simulations with saved model
        self.RL_mode = "test"
        self.model_name = "v00"

        # "ABOD": Only match on the major antigens.
        # "relimm": Use relative immunogenicity weights for mismatching.
        # "patgroups": Use patient group specific mismatching weights.
        self.strategy = "relimm"

        self.demand_scenario_names = [f"regional_{i}" for i in range(5)]
        self.supply_scenario_names = [f"cau{round(self.donor_eth_distr[0]*100)}_afr{round(self.donor_eth_distr[1]*100)}_asi{round(self.donor_eth_distr[2]*100)}" for i in range(5)]


        ########################################
        # ONLY USED FOR REINFORCEMENT LEARNING #

        # self.nn_update_iter = 5
        self.max_memory_size = 1000
        self.batch_size = 20
        self.gamma = 0.99
        self.alpha = 0.01
        self.dropout = 0.2
        self.exploration_max = 1.0
        self.exploration_min = 0.01
        self.exploration_decay = 0.1


        ###########################################
        # SETTINGS CONSIDERING BLOOD AND MATCHING #

        self.max_age = 35

        self.ABOD = ["O-", "O+", "A-", "A+", "B-", "B+", "AB-", "AB+"]
        self.major = ["A", "B", "D"]
        self.minor = [] # "C", "c", "E", "e", "K", "k", "M", "N", "S", "s", "Fya", "Fyb", "Jka", "Jkb"
        self.relimm = [0, 0, 0, 0.0345, 0.0705, 0.2395, 0.0836, 0.3838, 0, 0.0296, 0, 0.0131, 0, 0.0443, 0.0131, 0.0836, 0.0033]
        self.donor_ABOD_distr = {"O-":0.1551, "O+":0.3731, "A-":0.0700, "A+":0.3074, "B-":0.0158, "B+":0.0604, "AB-":0.0047, "AB+":0.0134}
        
        #A,B,D,C,c,E,e,K,k,M,N,S,s,Fya,Fyb,Jka,Jkb
        self.patgroup_weights = [
            ["M","M","M",1,1,1,1,1,1,0,0,1,1,1,1,1,1],                  # Other
            ["M","M","M",1,"M","M",1,"M",1,0,0,1,1,1,1,1,1],            # Wu45
            ["M","M","M","M","M","M","M","M",1,0,0,1,1,2,2,2,2],        # MDS
            ["M","M","M","M","M","M","M","M",1,0,0,1,1,2,2,2,2],        # Thal
            ["M","M","M","M","M","M","M","M",1,0,0,1,2,3,3,3,3],        # AIHA
            ["M","M","M","M","M","M","M","M",1,0,0,2,1,3,3,3,3],        # ALA
            ["M","M","M","M","M","M","M","M",1,0,0,3,2,"M",3,"M","M"],  # SCD
        ]

        self.request_lead_time_probabilities = {
            "Other" : [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 0, 0, 0, 0, 0, 0, 0],     # Other = uniform between [0,7] (CHANGE 3)
            "Wu45" : [1/2, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],                # Wu45 = 50/50 on the day or one day ahead.
            "MDS" : [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],                     # MDS = 7 days ahead
            "Thal" : [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],                    # Thal = 7 days ahead
            "AIHA" : [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],                    # AIHA = on the same day (CHANGE 5)
            "ALA" : [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 0, 0, 0, 0, 0, 0, 0],       # ALA = uniform between [0,7] (CHANGE 4)
            "SCD" : [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]                      # SCD = 7 days ahead
        }

        # each row specifies the probability of the corresponding patient type having a demand for [1,2,3,4] units respectively
        self.request_num_units_probabilities = {
            "Other" : [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782], 
            "Wu45" : [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782], 
            "MDS" : [0, 1, 0, 0], 
            "Thal" : [0, 1, 0, 0], 
            "AIHA" : [0, 1, 0, 0], 
            "ALA" : [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782], 
            "SCD" : [0, 1, 0, 0],
        }

        self.ABO_genotypes = [
            [ 0, 0 ],
            [ 0, 1 ],
            [ 1, 0 ],
            [ 1, 1 ]]

        self.ABO_prevalences = {
            "Caucasian" : [0.43, 0.09, 0.44, 0.04],
            "African" : [0.27, 0.49, 0.2 , 0.04],
            "Asian" : [0.27, 0.43, 0.25, 0.05]}

        self.Rhesus_genotypes = [
            [0,1,1,1,1],
            [0,1,0,1,0],
            [0,1,0,0,1],
            [0,0,1,1,0],
            [0,0,1,0,1],
            [0,1,1,1,0],
            [0,1,1,0,1],
            [0,1,0,1,1],
            [0,0,1,1,1],
            [1,1,1,1,1],
            [1,1,0,1,0],
            [1,1,0,0,1],
            [1,0,1,1,0],
            [1,0,1,0,1],
            [1,1,1,1,0],
            [1,1,1,0,1],
            [1,1,0,1,1],
            [1,0,1,1,1]]

        self.Rhesus_prevalences = {
            "Caucasian" : [0.   , 0.   , 0.   , 0.   , 0.151, 0.   , 0.008, 0.   , 0.009, 0.133, 0.   , 0.185, 0.023, 0.021, 0.001, 0.349, 0.002, 0.118],
            "African" : [0.   , 0.   , 0.   , 0.   , 0.068, 0.   , 0.   , 0.   , 0.   , 0.056, 0.   , 0.02 , 0.002, 0.458, 0.   , 0.21 , 0.   , 0.186],
            "Asian" : [0.   , 0.   , 0.001, 0.001, 0.001, 0.   , 0.001, 0.   , 0.   , 0.303, 0.   , 0.518, 0.044, 0.003, 0.004, 0.085, 0.014, 0.025]}

        self.Kell_genotypes = [
            [ 0, 0 ],
            [ 1, 0 ],
            [ 0, 1 ],
            [ 1, 1 ]]

        self.Kell_prevalences = {
            "Caucasian" : [0.   , 0.002, 0.91 , 0.088],
            "African" : [0.   , 0.   , 0.98 , 0.02 ],
            "Asian" : [0.   , 0.   , 1.   , 0.   ]}

        self.MNS_genotypes = [
            [1,0,1,0],
            [1,0,1,1],
            [1,0,0,1],
            [1,1,1,0],
            [1,1,1,1],
            [1,1,0,1],
            [0,1,1,0],
            [0,1,1,1],
            [0,1,0,1],
            [1,0,0,0],
            [1,1,0,0],
            [0,1,0,0]]

        self.MNS_prevalences  = {
            "Caucasian" : [0.06 , 0.14 , 0.08 , 0.04 , 0.24 , 0.22 , 0.01 , 0.06 , 0.15 , 0.   , 0.   , 0.   ],
            "African" : [0.02 , 0.07 , 0.16 , 0.02 , 0.13 , 0.325, 0.02 , 0.05 , 0.19 , 0.004, 0.004, 0.007],
            "Asian" : [0.06 , 0.14 , 0.08 , 0.04 , 0.24 , 0.22 , 0.01 , 0.06 , 0.15 , 0.   , 0.   , 0.   ]}

        self.Duffy_genotypes = [
            [ 0, 0 ],
            [ 1, 0 ],
            [ 0, 1 ],
            [ 1, 1 ]]

        self.Duffy_prevalences = {
            "Caucasian" : [0.   , 0.17 , 0.34 , 0.49 ],
            "African" : [0.68 , 0.09 , 0.22 , 0.01 ],
            "Asian" : [0.   , 0.908, 0.003, 0.089]}

        self.Kidd_genotypes = [
            [ 0, 0 ],
            [ 1, 0 ],
            [ 0, 1 ],
            [ 1, 1 ]]

        self.Kidd_prevalences = {
            "Caucasian" : [0.   , 0.263, 0.234, 0.503],
            "African" : [0.   , 0.511, 0.081, 0.488],
            "Asian" : [0.009, 0.232, 0.268, 0.491]}



#         PATIENT_DEMAND_PROBABILITIES = [
#     [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782],
#     [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782],
#     [0,1,0,0],
#     [0,1,0,0],
#     [0,1,0,0],
#     [0.40437368293541415, 0.4968390449938975, 0.06828851313979055, 0.03049875893089782],
#     [0,1,0,0],
# ]

# #specifies the probability per patient category for how many days in advance the request becomes kwown [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
# PATIENT_REQUEST_AVAILABLE_PROBABILITIES = [
#     [0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0],
#     [0.5,0.5,0,0,0,0,0,0,0,0,0,0,0,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0,0,0,0],
#     [1/7,1/7,1/7,1/7,1/7,1/7,1/7,0,0,0,0,0,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0,0,0,0]
# ]

# #specifies the probability of a patient group being a certain ethnicity (cau, neg, asi)
# PATIENT_TYPE_ETHNICITY_DISTRIBUTIONS = [
#     [1,0,0],
#     [1,0,0],
#     [1,0,0],
#     [1,0,0],
#     [1,0,0],
#     [1,0,0],
#     [0,1,0],
# ]




# #probabilities for getting a patient from a certain category
# OLVG = [
#     .886867,
#     .049352,
#     0,
#     .008607969,
#     .014933,
#     .031632,
#     .008607969
# ]

# AMC = [
#     .64605,
#     .10250,
#     .04542,
#     .05665,
#     .02731,
#     .06543,
#     .05665
# ]
# Normal = [
#     1,
#     0,
#     0,
#     0,
#     0,
#     0,
#     0
# ]