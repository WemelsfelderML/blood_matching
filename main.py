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

    # Sample demand for each day in the simulation, and write to a csv file.
    if SETTINGS.mode == "demand":
        for htype in SETTINGS.n_hospitals.keys():
            for _ in range(SETTINGS.n_hospitals[htype]):
                generate_demand(SETTINGS, PARAMS, htype, SETTINGS.avg_daily_demand[htype])

    # Sample RBC units to be used as supply in the simulation, and write to csv file.
    elif SETTINGS.mode == "supply":
        generate_supply(SETTINGS, PARAMS)

    # Run the simulation, using either linear programming or reinforcement learning to determine the matching strategy.
    elif SETTINGS.mode == "optimize":

        print(f"Starting {SETTINGS.method} optimization ({SETTINGS.line}line).")
        print(f"Results will be written to {SETTINGS.model_name} folder.")
        print(f"Simulating {SETTINGS.init_days + SETTINGS.test_days} days ({SETTINGS.init_days} init, {SETTINGS.test_days} test).")
        if sum(SETTINGS.n_hospitals.values()) > 1:
            print(f"Multi-hospital scenario with {SETTINGS.n_hospitals['regional']} regional and {SETTINGS.n_hospitals['university']} university hospitals.")
        else:
            print(f"Single-hospital scenario, {max(SETTINGS.n_hospitals, key = lambda i: SETTINGS.n_hospitals[i])} hospital.")
        print(f"Using {SETTINGS.strategy} strategy for matching, with patgroup_musts = {SETTINGS.patgroup_musts}.")

        if SETTINGS.method == "LP":
            simulation(SETTINGS, PARAMS)
        # elif SETTINGS.method == "RL":
        #     reinforcement_learning(SETTINGS, PARAMS)
        else:
            print("Parameter 'mode' is set to 'optimize', but no existing method for optimization is given. Try 'RL' or 'LP'.")
    else:
        print("No mode for running the code is given. Please change the 'mode' parameter in 'settings.py' to one of the following values:")
        print("'demand': generate demand scenarios")
        print("'supply': generate supply scenarios")
        print("'optimize': for optimizing RBC matching, either using LP or RL method")


# Check whether a given path exists, and create the path if it doesn't.
def check_dir_existence(path):
    if os.path.exists(path) == False:
        os.mkdir(path)


if __name__ == "__main__":

    main()