import pandas as pd
import os
import pickle

def save_state(SETTINGS, df, e, day, dc, hospitals):

	path = SETTINGS.home_dir + "wip"
	SETTINGS.check_dir_existence(path)
	path += f"/{SETTINGS.model_name}"
	SETTINGS.check_dir_existence(path)
	path += f"/{e}"
	SETTINGS.check_dir_existence(path)
	path += f"/{SETTINGS.strategy}_{'-'.join([str(SETTINGS.n_hospitals[ds]) + ds[:3] for ds in SETTINGS.n_hospitals.keys()])}"	


	df.to_csv(path + "_df.csv", sep=',', index=True)
	dc.pickle(path + "_dc.pickle")
	for h in range(len(hospitals)):
		hospitals[h].pickle(path + f"_h{h}.pickle")


def load_state(SETTINGS, e, df, dc, hospitals):

	path = SETTINGS.home_dir + f"wip/{SETTINGS.model_name}/{e}/{SETTINGS.strategy}_{'-'.join([str(SETTINGS.n_hospitals[ds]) + ds[:3] for ds in SETTINGS.n_hospitals.keys()])}"	

	if os.path.exists(path + "_df.csv") == True:

		df = pd.read_csv(path + "_df.csv")
		day = max(df[df["logged"]==True]["day"]) + 1
		df = df.set_index(["day", "location"])
		
		dc = unpickle(path + "_dc.pickle")

		hospitals = []
		for h in range(sum(SETTINGS.n_hospitals.values())):
			hospitals.append(unpickle(path + f"_h{h}.pickle"))

	else:
		day = 0

	return df, day, dc, hospitals


def unpickle(path):
	with open(path, 'rb') as f:
		return pickle.load(f)