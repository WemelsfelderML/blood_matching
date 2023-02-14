import numpy as np
import re

# Take a solved MINRAR model and get the value for each of the variables.
def read_minrar_solution(SETTINGS, PARAMS, df, model, dc, hospitals, episode, day=0):

    # Variables for the multi-hospital scenario.
    if sum(SETTINGS.n_hospitals.values()) > 1:

        # Create numpy arrays filled with zeros.
        xh = [np.zeros([len(hospital.inventory), len(hospital.requests)]) for hospital in hospitals]
        xdc = [np.zeros([len(dc.inventory), len(hospital.requests)]) for hospital in hospitals]
        y = [np.zeros([len(hospital.requests)]) for hospital in hospitals]
        z = [np.zeros([len(hospital.requests), len(PARAMS.major + PARAMS.minor)]) for hospital in hospitals]

        # Get the values for all model variables and write them into the numpy arrays created above.
        for h in range(len(hospitals)):
            for var in model.getVars():
                var_name = re.split(r'\W+', var.varName)[0]
                if var_name == f"xh{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xh[h][index0, index1] = var.X
                if var_name == f"xdc{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    xdc[h][index0, index1] = var.X
                if var_name == f"y{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    y[h][index0] = var.X
                if var_name == f"z{h}":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    z[h][index0, index1] = var.X

        # Calculate the number of variables and add this information to the output dataframe.
        nvars = sum([sum([np.product(var.shape) for var in [xh[h], xdc[h], y[h], z[h]]]) for h in range(len(hospitals))])
        print("nvars:",nvars)
        df.loc[(day,hospitals[h].name),"nvars"] = nvars

        return df, xh, xdc, y, z

    # Variables for the single-hospital scenario.
    else:

        # Create numpy arrays filled with zeros.
        x = np.zeros([len(hospitals[0].inventory), len(hospitals[0].requests)])
        y = np.zeros([len(hospitals[0].requests)])
        z = np.zeros([len(hospitals[0].requests), len(PARAMS.major + PARAMS.minor)])

        # Get the values for all model variables and write them into the numpy arrays created above.
        for var in model.getVars():
            var_name = re.split(r'\W+', var.varName)[0]
            if var_name == "x":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                x[index0, index1] = var.X
            if var_name == "y":
                index0 = int(re.split(r'\W+', var.varName)[1])
                y[index0] = var.X
            if var_name == "z":
                index0 = int(re.split(r'\W+', var.varName)[1])
                index1 = int(re.split(r'\W+', var.varName)[2])
                z[index0, index1] = var.X

        # Do the same for two other variables in case the offline optimization was performed.
        if SETTINGS.line == "off":
            a = np.zeros([len(hospitals[0].inventory), SETTINGS.init_days + SETTINGS.test_days])
            b = np.zeros([len(hospitals[0].inventory), SETTINGS.init_days + SETTINGS.test_days])
            for var in model.getVars():
                var_name = re.split(r'\W+', var.varName)[0]
                if var_name == "a":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    a[index0, index1] = var.X
                if var_name == "b":
                    index0 = int(re.split(r'\W+', var.varName)[1])
                    index1 = int(re.split(r'\W+', var.varName)[2])
                    b[index0, index1] = var.X

            # Calculate the number of variables and add this information to the output dataframe.
            nvars = sum([np.product(var.shape) for var in [x, y, z, a, b]])
            print("nvars:",nvars)
            df.loc[(day,hospitals[0].name),"nvars"] = nvars

            return df, x, y, z, a, b

        else:
            # Calculate the number of variables and add this information to the output dataframe.
            nvars = max(df.loc[(day,hospitals[0].name),"nvars"], sum([np.product(var.shape) for var in [x, y, z]]))
            print("nvars:",nvars)
            df.loc[(day,hospitals[0].name),"nvars"] = nvars

            return df, x, y, z

    # Write the values found to pickle files.
    # with open(SETTINGS.generate_filename("results") + f"x_{SETTINGS.strategy}_{hospital.htype[:3]}_{episode}-{day}.pickle", "wb") as f:
    #     pickle.dump(x, f)
    # with open(SETTINGS.generate_filename("results") + f"y_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(y, f)
    # with open(SETTINGS.generate_filename("results") + f"z_{SETTINGS.strategy}_{hospital.htype[:3]}_{e}.pickle", "wb") as f:
    #     pickle.dump(z, f)


# Take a solved MINRAR model and abstract the optimal variable values.
def read_transports_solution(model, size_I, size_H):

    # Get the values of the model variable as found for the optimal solution.
    x = np.zeros([size_I, size_H])
    for var in model.getVars():
        var_name = re.split(r'\W+', var.varName)[0]
        if var_name == "x":
            index0 = int(re.split(r'\W+', var.varName)[1])
            index1 = int(re.split(r'\W+', var.varName)[2])
            x[index0, index1] = var.X

    return x

