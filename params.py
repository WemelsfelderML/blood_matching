class Params():
    
    def __init__(self, SETTINGS):

        ####################
        # BLOOD PARAMETERS #
        ####################

        self.max_age = 35

        self.ABOD = ["O-", "O+", "A-", "A+", "B-", "B+", "AB-", "AB+"]
        self.major = ["A", "B", "D"]
        self.minor = ["C", "c", "E", "e", "K", "k", "M", "N", "S", "s", "Fya", "Fyb", "Jka", "Jkb"] # "C", "c", "E", "e", "K", "k", "M", "N", "S", "s", "Fya", "Fyb", "Jka", "Jkb"
        
        self.relimm_weights = [0, 0, 0, 0.0345, 0.0705, 0.2395, 0.0836, 0.3838, 0, 0.0296, 0, 0.0131, 0, 0.0443, 0.0131, 0.0836, 0.0033]


        ##################
        # PATIENT GROUPS #
        ##################

        self.patgroups = ["Other", "Wu45", "MDS", "Thal", "AIHA", "ALA", "SCD"]

        if SETTINGS.demand_scenario == "Other":
            self.patgroup_distr = [1, 0, 0, 0, 0, 0, 0]
        elif SETTINGS.demand_scenario == "regional":
            self.patgroup_distr = [0.886867, 0.049352, 0, 0.008607969, 0.014933, 0.031632, 0.008607969]
        elif SETTINGS.demand_scenario == "university":
            self.patgroup_distr = [0.64605, 0.10250, 0.04542, 0.05665, 0.02731, 0.06543, 0.05665]
        
        # A,B,D,C,c,E,e,K,k,M,N,S,s,Fya,Fyb,Jka,Jkb
        # self.patgroup_weights = [
        #     ["M","M","M",1,1,1,1,1,1,0,0,1,1,1,1,1,1],                  # Other
        #     ["M","M","M",1,"M","M",1,"M",1,0,0,1,1,1,1,1,1],            # Wu45
        #     ["M","M","M","M","M","M","M","M",1,0,0,1,1,2,2,2,2],        # MDS
        #     ["M","M","M","M","M","M","M","M",1,0,0,1,1,2,2,2,2],        # Thal
        #     ["M","M","M","M","M","M","M","M",1,0,0,1,2,3,3,3,3],        # AIHA
        #     ["M","M","M","M","M","M","M","M",1,0,0,2,1,3,3,3,3],        # ALA
        #     ["M","M","M","M","M","M","M","M",1,0,0,3,2,"M",3,"M","M"],  # SCD
        # ]
        self.patgroup_weights = [
            [-1, -1, -1, 0.02651, 0.05429, 0.18433, 0.06439, 0.29543, 0, 0, 0, 0.00337, 0.00128, 0.02273, 0.00673, 0.04293, 0.00168],   # Other
            [-1, -1, -1, 0.02651, -1, -1, 0.06439, -1, 0, 0, 0, 0.00337, 0.00128, 0.02273, 0.00673, 0.04293, 0.00168],                  # Wu45
            [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0.10100, 0.03844, 0.68177, 0.20200, 1.28778, 0.05050],                            # MDS
            [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0.10100, 0.03844, 0.68177, 0.20200, 1.28778, 0.05050],                            # Thal
            [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0.20200, 0.07689, 01.3653, 0.40401, 2.57556, 0.10100],                            # AIHA
            [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0.20200, 0.07689, 01.3653, 0.40401, 2.57556, 0.10100],                            # ALA
            [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0.33667, 0.12815, -1, 0.67335, -1, -1],                                           # SCD
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

        self.donor_ABOD_distr = {"O-":0.1551, "O+":0.3731, "A-":0.0700, "A+":0.3074, "B-":0.0158, "B+":0.0604, "AB-":0.0047, "AB+":0.0134}

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

        
        ##########
        # GUROBI #
        ##########
        
        self.status_code = {
            1 : "MODEL IS LOADED, BUT NO SOLUTION IS AVAILABLE",
            2 : "MODEL IS SOLVED OPTIMALLY",
            3 : "MODEL IS INFEASIBLE",
            4 : "MODEL IS EITHER INFEASIBLE OR UNBOUNDED\nTo obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.",
            5 : "MODEL IS UNBOUNDED\nAn unbounded ray allows the objective to improve without limit.",
            6 : "NO SOLUTION AVAILABLE\nThe optimal objective was worse than the Cutoff parameter.",
            7 : "OPTIMIZATION TERMINATED\nIterationLimit or BarIterLimit parameter was exceeded.",
            8 : "OPTIMIZATION TERMINATED\nNodeLimit parameter was exceeded.",
            9 : "OPTIMIZATION TERMINATED\nTimeLimit parameter was exceeded.",
            10 : "OPTIMIZATION TERMINATED\nSolutionLimit parameter was exceeded.",
            11 : "OPTIMIZATION TERMINATED BY USER\nObtained results are not saved.",
            12 : "OPTIMIZATION TERMINATED\nUnrecoverable numerical difficulties.",
            13 : "UNABLE TO SATISFY OPTIMALITY\nA sub-optimal solution is available.",
            14 : "ASYNCHRONOUS CALL WAS MADE, ASSOCIATED OPTIMIZATION NOT YET COMPLETE",
            15 : "LIMIT SET BY USER WAS EXCEEDED\nThis is either a bound on the best objective or the best bound."
        }



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