import random
import numpy as np
from operator import itemgetter

# Class representing the phenotype of an individual
class Blood:

    def __init__(self, PARAMS, ethnicity = None, patgroup = "", major = "", minor = "", num_units = 0, issuing_day = 0, age = 0):

        vector = []

        if major == "":
            vector += random.choices(PARAMS.ABO_genotypes, weights = PARAMS.ABO_prevalences[ethnicity], k=1)[0]
            vector += random.choices(PARAMS.Rhesus_genotypes, weights = PARAMS.Rhesus_prevalences[ethnicity], k=1)[0]

        else:

            if "A" in major:
                vector += [1]
            else:
                vector += [0]
            
            if "B" in major:
                vector += [1]
            else:
                vector += [0]

            if minor == "":
                Rhesus_genotypes_Dpos = []
                Rhesus_genotypes_Dneg = []
                Rhesus_prevalences_Dpos = []
                Rhesus_prevalences_Dneg = []

                for i in range(len(PARAMS.Rhesus_genotypes)):
                    genotype = PARAMS.Rhesus_genotypes[i]
                    prevalence = PARAMS.Rhesus_prevalences[ethnicity][i]
                    if genotype[0] == 1:
                        Rhesus_genotypes_Dpos.append(genotype)
                        Rhesus_prevalences_Dpos.append(prevalence)
                    else:
                        Rhesus_genotypes_Dneg.append(genotype)
                        Rhesus_prevalences_Dneg.append(prevalence)

                if "+" in major:
                    vector += random.choices(Rhesus_genotypes_Dpos, weights = Rhesus_prevalences_Dpos, k=1)[0]
                elif "-" in major:
                    vector += random.choices(Rhesus_genotypes_Dneg, weights = Rhesus_prevalences_Dneg, k=1)[0]
                else:
                    print("Attempted to create Blood instance, but the given major bloodgroup does not contain + or - for RhD.")

            else:
                if "+" in major:
                    vector += [1]
                elif "-" in major:
                    vector += [0]
                else:
                    print("Attempted to create Blood instance, but the given major bloodgroup does not contain + or - for RhD.")

        if minor == "":
            vector += random.choices(PARAMS.Kell_genotypes, weights = PARAMS.Kell_prevalences[ethnicity], k=1)[0]
            vector += random.choices(PARAMS.MNS_genotypes, weights = PARAMS.MNS_prevalences[ethnicity], k=1)[0]
            vector += random.choices(PARAMS.Duffy_genotypes, weights = PARAMS.Duffy_prevalences[ethnicity], k=1)[0]
            vector += random.choices(PARAMS.Kidd_genotypes, weights = PARAMS.Kidd_prevalences[ethnicity], k=1)[0]

            antigens_vector = ["A", "B", "D", "C", "c", "E", "e", "K", "k", "M", "N", "S", "s", "Fya", "Fyb", "Jka", "Jkb"]
            self.vector = [vector[antigens_vector.index(ag)] for ag in (PARAMS.major + PARAMS.minor)]

        else:
            vector += minor
            self.vector = vector

        self.major = vector_to_major(vector)

        # Only used for inventory products.
        self.ethnicity = ethnicity
        self.age = age

        # Only used for requests.
        self.patgroup = patgroup
        self.num_units = num_units
        self.issuing_day = issuing_day
        self.best_mismatch_penalty = 999
        self.allocated_from_dc = 0
        

    # Transform the binary antigen vector to a blood group index
    def vector_to_bloodgroup_index(self):
        return int("".join(str(i) for i in self.vector),2)


    def get_usability(self, PARAMS, hospitals, antigens = []):

        # IMPORTANT: this is now hardcoded for the case where SCD patients are Africans and all others are Caucasions.
        avg_daily_demand_african = sum([PARAMS.patgroup_distr[hospital.htype]["SCD"] * hospital.avg_daily_demand for hospital in hospitals])
        avg_daily_demand_total = sum([hospital.avg_daily_demand for hospital in hospitals])
        part_african = avg_daily_demand_african / avg_daily_demand_total

        if antigens == []:
            usability_ABO = 0
            usability_RhD = 1

            # Calculate the ABO-usability of this blood product, by summing all prevalences of the phenotypes that can receive this product.
            ABO_v = self.vector[:2] 
            ABO_g = PARAMS.ABO_genotypes
            for g in range(len(ABO_g)):
                if all(v <= g for v, g in zip(ABO_v, ABO_g[g])):
                    usability_ABO += PARAMS.ABO_prevalences["African"][g] * part_african
                    usability_ABO += PARAMS.ABO_prevalences["Caucasian"][g] * (1 - part_african)

            # Calculate the RhD-usability of this blood product, by summing all prevalences of the phenotypes that can receive this product.
            # If the considered blood product is RhD negative, usability is always 1. Therefore usability is only calculated when the product is RhD positive.
            # TODO make more efficient by pre-computing these values and storing them somewhere
            Dpos = np.array([g[0] for g in PARAMS.Rhesus_genotypes])
            Dpos_prevalence = sum(np.array(PARAMS.Rhesus_prevalences["African"]) * part_african * Dpos) + sum(np.array(PARAMS.Rhesus_prevalences["Caucasian"]) * (1 - part_african) * Dpos)
            if self.vector[2] == 1:
                usability_RhD = Dpos_prevalence

            #return the product of all the invdiviual system usabilities to compute the final usabilty
            return usability_ABO * usability_RhD

        else:
            antigens = [ag for ag in (PARAMS.major + PARAMS.minor) if ag in antigens]   # Get intersection of all antigens given to consider, and all antigens in the model.
            usability_ABO = self.get_usability_system(["A", "B"], antigens, PARAMS.ABO_genotypes, PARAMS.ABO_prevalences, part_african)
            usability_Rhesus = self.get_usability_system(["D", "C", "c", "E", "e"], antigens, PARAMS.Rhesus_genotypes, PARAMS.Rhesus_prevalences, part_african)
            usability_Kell = self.get_usability_system(["K", "k"], antigens, PARAMS.Kell_genotypes, PARAMS.Kell_prevalences, part_african)
            usability_MNS = self.get_usability_system(["M", "N", "S", "s"], antigens, PARAMS.MNS_genotypes, PARAMS.MNS_prevalences, part_african)
            usability_Duffy = self.get_usability_system(["Fya", "Fyb"], antigens, PARAMS.Duffy_genotypes, PARAMS.Duffy_prevalences, part_african)
            usability_Kidd = self.get_usability_system(["Jka", "Jkb"], antigens, PARAMS.Kidd_genotypes, PARAMS.Kidd_prevalences, part_african)

        #return the product of all the invdiviual system usabilities to compute the final usabilty
        return usability_ABO * usability_Rhesus * usability_Kell * usability_MNS * usability_Duffy * usability_Kidd


    def get_usability_system(self, system_antigens, antigens, genotypes, prevalences, part_african):

        # TODO: now the usability is only calculated if all antigens of the system are included. Extend this to calculating it for only selected antigens.
        if all(ag in antigens for ag in system_antigens):
            
            usability = 0
            vector_indices = [antigens.index(k) for k in system_antigens]

            # Calculate the ABO-usability of this blood product, by summing all prevalences of the phenotypes that can receive this product.
            vector = [self.vector[i] for i in vector_indices]
            for g in range(len(genotypes)):
                if all(v <= g for v, g in zip(vector, genotypes[g])):
                    usability += prevalences["African"][g] * part_african
                    usability += prevalences["Caucasian"][g] * (1 - part_african)

            return usability

        else:
            return 1
            



# Obtain the major blood group from a blood antigen vector
def vector_to_major(vector):

    major = ""

    if vector[0] == 1:
        major += "A"
    if vector[1] == 1:
        major += "B"
    if len(major) == 0:
        major += "O"
    if vector[2] == 1:
        major += "+"
    else:
        major += "-"

    return major

def timewise_possible(SETTINGS, PARAMS, I, R, day):
    
    T = np.zeros([len(I), len(R)])
    for i in I.keys():
        for r in R.keys():
            T[i,r] = 1 if (PARAMS.max_age - 1 - I[i].age) >= (R[r].issuing_day - day) else 0

    return T

def precompute_compatibility(SETTINGS, PARAMS, I, R):

    antigens = PARAMS.major + PARAMS.minor
    C = np.zeros([len(I), len(R)])

    if ("patgroups" in SETTINGS.strategy) or SETTINGS.patgroup_musts:
        for i in range(len(I)):
            for r in range(len(R)):
                v_musts_ir = [(I[i].vector[k], R[r].vector[k]) for k in range(len(antigens)) if PARAMS.patgroup_weights.loc[R[r].patgroup,antigens[k]] == 10]
                if all(vi <= vr for vi, vr in v_musts_ir):
                    C[i,r] = 1

    else:
        num_major = len(PARAMS.major)
        for i in range(len(I)):
            for r in range(len(R)):
                if all(vi <= vr for vi, vr in zip(I[i].vector[:num_major], R[r].vector[:num_major])):
                    C[i,r] = 1

    return C

def compatibility(SETTINGS, PARAMS, ip, rq):

    antigens = PARAMS.major + PARAMS.minor

    if ("patgroups" in SETTINGS.strategy) or SETTINGS.patgroup_musts:
        v_musts_ir = [(ip.vector[k], rq.vector[k]) for k in range(len(antigens)) if PARAMS.patgroup_weights.loc[rq.patgroup,antigens[k]] == 10]
        if all(vi <= vr for vi, vr in v_musts_ir):
            return 1
        else:
            return 0

    else:
        num_major = len(PARAMS.major)
        if all(vi <= vr for vi, vr in zip(ip.vector[:num_major], rq.vector[:num_major])):
            return 1
        else:
            return 0



#         COMPATIBILITY_DONOR_PATIENT = [
#             [1,1,1,1,1,1,1,1],
#             [0,1,0,1,0,1,0,1],
#             [0,0,1,1,0,0,1,1],
#             [0,0,0,1,0,0,0,1],
#             [0,0,0,0,1,1,1,1],
#             [0,0,0,0,0,1,0,1],
#             [0,0,0,0,0,0,1,1],
#             [0,0,0,0,0,0,0,1]]

#         COMPATIBILITY_PATIENT_DONOR = [
#             [1,0,0,0,0,0,0,0],
#             [1,1,0,0,0,0,0,0],
#             [1,0,1,0,0,0,0,0],
#             [1,1,1,1,0,0,0,0],
#             [1,0,0,0,1,0,0,0],
#             [1,1,0,0,1,1,0,0],
#             [1,0,1,0,1,0,1,0],
#             [1,1,1,1,1,1,1,1]]



#         #/ <summary>
#         #/ Class representing a blood type of an individual
#         #/ </summary>
#         #/ <param name="vector">Binary integer array indicating the presence of the antigens</param>
#         #/ <param name="ethnicity"></param>
# #C# TO PYTHON CONVERTER TODO TASK: There is no Python equivalent to multiple constructors:
# #ORIGINAL LINE: public Blood(int[] vector, Ethnicity ethnicity)
#         def __init__(self, vector, ethnicity):
#             self._initialize_instance_fields()
#             self.Ethnicity = ethnicity
#             self.Vector = vector
#             self.Major = VectorToMajor(vector)

#         

    # import math

# #/ <summary>
# #/ Precompute the different usability values based on the ethnic distribution of patient in the given demand scenario
# #/ </summary>
# #/ <param name="demand"></param>
# def PrecomputeUsabilities(self, demand):
#     #precompute different usabilities
#     self.Usability_a3 = BloodGroupInfo.GetUsability(A3, Vector, demand)
#     self.Usability_a8 = BloodGroupInfo.GetUsability(A8, Vector, demand)
#     self.Usability_a14 = BloodGroupInfo.GetUsability(A14, Vector, demand)
#     self.Usability_aminor = BloodGroupInfo.GetUsability(AMinor, Vector, demand)
#     usabilityIsPrecomputed = True

# def PrecomputeUsabilitiesDC(self, demandScenarios):
#     #precompute different usabilities
#     Ua3 = 0.0
#     Ua8 = 0.0
#     Ua14 = 0.0
#     Uaminor = 0.0
#     for demand in demandScenarios:
#         Ua3 += BloodGroupInfo.GetUsability(A3, Vector, demand)
#         Ua8 += BloodGroupInfo.GetUsability(A8, Vector, demand)
#         Ua14 += BloodGroupInfo.GetUsability(A14, Vector, demand)
#         Uaminor += BloodGroupInfo.GetUsability(AMinor, Vector, demand)

#     self.Usability_a3 = Ua3 / len(demandScenarios)
#     self.Usability_a8 = Ua8 / len(demandScenarios)
#     self.Usability_a14 = Ua14 / len(demandScenarios)
#     self.Usability_aminor = Uaminor / len(demandScenarios)

#     usabilityIsPrecomputed = True

# #Usability of this blood with respect to a certain antigenset
# def Usability(self, antigen_set):
#     if not usabilityIsPrecomputed:
#         raise Exception("Usability values have not yet been precomputed! Please call PrecomputeUsabilities(demandScenario) before accessing these variables")


#     if antigen_set is AntigenSet.A3:
#         return self.Usability_a3
#     if antigen_set is AntigenSet.A8:
#         return self.Usability_a8
#     if antigen_set is AntigenSet.A14:
#         return self.Usability_a14
#     if antigen_set is AntigenSet.AMinor:
#         return self.Usability_aminor
#     if antigen_set is AntigenSet.None_:
#         return 0
#     return -1



# #/ <summary>
# #/ Return major blood group compatibility as binary integer
# #/ </summary>
# #/ <param name="from"></param>
# #/ <param name="to"></param>
# #/ <returns></returns>
# @staticmethod
# def MajorCompatible(from_, to):
#     if math.fmod(int(from_), 2) >math.fmod(int(to), 2):
#         return 0
#     if int(from_) > int(to):
#         return 0
#     if (from_ is Major.Aneg or from_ is Major.Apos) and (to is Major.Bneg or to is Major.Bpos):
#         return 0
#     return 1


# #/ <summary>
# #/ Get compatibility between two blood type vectors based on a set of antigens
# #/ </summary>
# #/ <param name="vec_from"></param>
# #/ <param name="vec_to"></param>
# #/ <param name="antigen_set"></param>
# #/ <returns></returns>
# @staticmethod
# def IsCompatible(vec_from, vec_to, antigen_set):
#     comp = True
#     set = A3
#     if antigen_set is AntigenSet.A3:
#         set = A3
#     if antigen_set is AntigenSet.A8:
#         set = A8
#     if antigen_set is AntigenSet.A14:
#         set = A14

#     k = 0
#     while k < len(vec_from):
#         comp = comp and (vec_from[k] <= vec_to[k] or (not set[k]))
#         k += 1
#     return 1 if comp else 0

