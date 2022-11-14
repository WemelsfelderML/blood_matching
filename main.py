from settings import *
from demand import *
from supply import *
from optimize import *

def main():

    SETTINGS = Settings()
    
    if SETTINGS.mode == "demand":
        generate_demand(SETTINGS)

    elif SETTINGS.mode == "supply":
        generate_supply(SETTINGS)

    elif SETTINGS.mode == "optimize":
        if SETTINGS.method == "RL":
            reinforcement_learning(SETTINGS)
        elif SETTINGS.method == "ILP":
            minrar(SETTINGS)
        else:
            print("Parameter 'mode' is set to 'optimize', but no existing method for optimization is given. Try 'RL' or 'ILP'.")
    else:
        print("No mode for running the code is given. Please change the 'mode' parameter in 'settings.py' to one of the following values:")
        print("'demand': generate demand scenarios")
        print("'supply': generate supply scenarios")
        print("'optimize': for optimizing RBC matching, either using ILP or RL method")


if __name__ == "__main__":


    main()