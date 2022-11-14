    # class PatientType(Enum):
    #     OTHER = 0
    #     WU45 = 1
    #     MDS = 2
    #     THAL = 3
    #     AIHA = 4
    #     ALA = 5      
    #     SCD = 6

    # class PatientTypeDistributionOrigin(Enum):
    #     ALLNORMAL = 0
    #     OLVG = 1
    #     AMC = 2

    


        #variables used for multi-hospital setup

#         # -1 is must
#         PatientTypeCostMatrix = [[-1, -1, -1, .02651, .05429, .18433, .06439, .29543, 0, 0, 0, .00337, .00128, .02273, .00673, .04293, .00168], 
#         [-1, -1, -1, .02651, -1, -1, .06439, -1, 0, 0, 0, .00337, .00128, .02273, .00673, .04293, .00168], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .10100, .03844, .68177, .20200, 1.28778, .05050], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .10100, .03844, .68177, .20200, 1.28778, .05050], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .20200, .07689, 1.3653, .40401, 2.57556, .10100], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .20200, .07689, 1.3653, .40401, 2.57556, .10100], 
#         [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, .33667, .12815, -1, .67335, -1, -1]]

        

#         #each of these arrays specifies the probability that a patient is of a specific category (in order of the PatientType enum above)
#         PatientTypeProbabilitiesOLVG = [.886867, .049352, 0, .008607969, .014933, .031632, .008607969]
#         PatientTypeProbabilitiesAMC = [.64605, .10250, .04542, .05665, .02731, .06543, .05665]
#         PatientTypeProbabilitiesAllNormal = [1, 0, 0, 0, 0, 0, 0]

#         #default is allnormal
#         PatientTypeProbabilities = PatientTypeProbabilitiesAllNormal
#         patientTypeDistributionOrigin = PatientTypeDistributionOrigin.ALLNORMAL

#         #/ <summary>
        

#         #/ Specifies the probability per patient type of the different lead times [0,14]
        


#         #/ <summary>
#         #/ Check whether a unit is compatible with a request based on the PatientTypeCostMatrix above
#         #/ </summary>
#         #/ <param name="unit"></param>
#         #/ <param name="request"></param>
#         #/ <returns></returns>
#         @staticmethod
#         def IsCompatible(unit, request):
#             type_index = int(request.PatientType)
#             k = 0
#             while k < self.Blood.Antigens.Length:
#                 if RBCallocation.Request.PatientTypeCostMatrix[type_index][k] == -1 and unit.Blood.Vector[k] == 1 and request.Blood.Vector[k] == 0:
#                     return 0
#                 k += 1
#             return 1


#         #/ <summary>
#         #/ Patient's request for blood
#         #/ </summary>
#         #/ <param name="blood"></param>
#         #/ <param name="day_needed"></param>
#         #/ <param name="day_available"></param>
#         #/ <param name="num_units"></param>
#         #/ <param name="patient_type"></param>
#         def __init__(self, blood, day_needed, day_available, num_units, patient_type):
#             #instance fields found by C# to Python Converter:
#             self.DayNeeded = 0
#             self.DayAvailable = 0
#             self.num_units = 0
#             self.Reciprocal = 0
#             self.PatientType = 0
#             self.Blood = None
#             self.MaximumMismatchPenalty = 100
#             self.HospitalIndex = 0
#             self.AllocatedAtDist = False
#             self.Index = 0

#             self.Blood = blood
#             self.DayAvailable = day_available
#             self.DayNeeded = day_needed
#             self.num_units = num_units
#             self.Reciprocal = 1 / float(num_units)
#             self.PatientType = patient_type


#         @staticmethod
#         def FromCsv(csvLine):
#             splitted = csvLine.split(',').ToList()
#             day_needed = int.Parse(splitted[0])
#             day_available = int.Parse(splitted[1])
#             num_units = int.Parse(splitted[2])
#             patient_type = None
#             if splitted[3] == "Wu45":
#                 patient_type = PatientType.WU45
#             elif splitted[3] != "Other":
#                 patient_type = Enum.Parse(typeof(PatientType), splitted[3])
#             else:
#                 patient_type = Enum.Parse(typeof(PatientType), splitted[3])
#             ethnicity = Enum.Parse(typeof(Ethnicity), splitted[4])

#             vector = Array.ConvertAll(splitted.GetRange(5, 17).ToArray(), int.Parse)
#             blood = Blood(vector, ethnicity)
#             return Request(blood, day_needed, day_available, num_units, patient_type)

#         @staticmethod
#         def ToCSVString(r):
#             result = [None for _ in range(5 + 17)]
#             result[0] = str(r.DayNeeded)
#             result[1] = str(r.DayAvailable)
#             result[2] = str(r.NumUnits)
#             result[3] = str(r.PatientType)
#             result[4] = str(r.Blood.Ethnicity)
#             for i in range(0, 17):
#                 result[5 + i] = str(r.Blood.Vector[i])
#             return string.Join(",", result)


#         #/ <summary>
#         #/ Sample a random patient request
#         #/ </summary>
#         #/ <param name="Program.Random"></param>
#         #/ <param name="day_needed"></param>
#         #/ <returns></returns>
#         @staticmethod
        
# }
