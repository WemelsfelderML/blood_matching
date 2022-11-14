import random

# Class representing the phenotype of an individual
class Blood:

    def __init__(self, SETTINGS = None, ethnicity = None, major = "", minor = ""):

        vector = []

        if major == "":
            vector += random.choices(SETTINGS.ABO_genotypes, weights = SETTINGS.ABO_prevalences[ethnicity], k=1)[0]
            vector += random.choices(SETTINGS.Rhesus_genotypes, weights = SETTINGS.Rhesus_prevalences[ethnicity], k=1)[0]

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

                for i in range(len(SETTINGS.Rhesus_genotypes)):
                    genotype = SETTINGS.Rhesus_genotypes[i]
                    prevalence = SETTINGS.Rhesus_prevalences[ethnicity][i]
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
            vector += random.choices(SETTINGS.Kell_genotypes, weights = SETTINGS.Kell_prevalences[ethnicity], k=1)[0]
            vector += random.choices(SETTINGS.MNS_genotypes, weights = SETTINGS.MNS_prevalences[ethnicity], k=1)[0]
            vector += random.choices(SETTINGS.Duffy_genotypes, weights = SETTINGS.Duffy_prevalences[ethnicity], k=1)[0]
            vector += random.choices(SETTINGS.Kidd_genotypes, weights = SETTINGS.Kidd_prevalences[ethnicity], k=1)[0]

        else:
            vector += minor

        self.ethnicity = ethnicity
        self.vector = vector
        self.major = vector_to_major(vector)


    # Transform the binary antigen vector to a blood group index
    def vector_to_bloodgroup_index(self):
        return int("".join(str(i) for i in self.vector),2)


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

