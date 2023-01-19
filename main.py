from settings import *
from params import *
from demand import *
from supply import *
from simulation import *
# from reinforcement_learning import *

def main():

    SETTINGS = Settings()
    PARAMS = Params(SETTINGS)

    # If a directory to store log files or results does not yet exist, make one.
    for path in ["results", "results/"+SETTINGS.model_name]:
        check_dir_existence(path)

    if SETTINGS.mode == "demand":
        for htype in SETTINGS.n_hospitals.keys():
            for _ in range(SETTINGS.n_hospitals[htype]):
                generate_demand(SETTINGS, PARAMS, htype, SETTINGS.avg_daily_demand[htype])

    elif SETTINGS.mode == "supply":
        generate_supply(SETTINGS, PARAMS)

    elif SETTINGS.mode == "optimize":

        if SETTINGS.method == "ILP":
            simulation(SETTINGS, PARAMS)
        # elif SETTINGS.method == "RL":
        #     reinforcement_learning(SETTINGS, PARAMS)
        else:
            print("Parameter 'mode' is set to 'optimize', but no existing method for optimization is given. Try 'RL' or 'ILP'.")
    else:
        print("No mode for running the code is given. Please change the 'mode' parameter in 'settings.py' to one of the following values:")
        print("'demand': generate demand scenarios")
        print("'supply': generate supply scenarios")
        print("'optimize': for optimizing RBC matching, either using ILP or RL method")


def check_dir_existence(path):
    if os.path.exists(path) == False:
        os.mkdir(path)

if __name__ == "__main__":


    main()