"""
Load in the data.
"""

from load_synergies import Synergies
from split_data import StratifiedKFolds
import matplotlib.pyplot as plt, numpy

# EXAMPLE: Load in data, but only those with QA = 1
synergies = Synergies()
synergies.process_data(remove_low_quality=True)
synergies_list = synergies.list

# EXAMPLE: Load in data, but select only cell lines from a given list
cell_lines = ['BT-20','CAL-120']
synergies_2 = Synergies()
synergies_2.process_data(remove_low_quality=True,select_cell_lines=cell_lines)

# EXAMPLE: Load in data, but select only drug combinations from a given list
drug_combinations = ['ADAM17.MTOR_1','ADAM17.AKT']
synergies_3 = Synergies()
synergies_3.process_data(remove_low_quality=True,select_drug_combinations=drug_combinations)

# EXAMPLE: Split data into test and train folds
kfolds = StratifiedKFolds(data_synergies=synergies.raw_synergies)
kfolds.process_data(remove_low_quality=True)
kfolds.split(folds=5,seed=1)
train, test = kfolds.return_split(index=0)

# EXAMPLE: Split data into test and train folds, but remove cell lines from a list because they have a single entry
remove_cell_lines = ['22RV1','VCaP']
kfolds = StratifiedKFolds(data_synergies=synergies.raw_synergies)
kfolds.process_data(remove_low_quality=True,remove_cell_lines=remove_cell_lines)
kfolds.split(folds=5,seed=1)
train, test = kfolds.return_split(index=0)
for cell_line in remove_cell_lines:
    assert cell_line not in kfolds.synergies.unique_cell_lines, "Removing cell lines not working!"

# EXAMPLE: split data and run cross-validation
kfolds_2 = StratifiedKFolds(data_synergies=synergies.raw_synergies)
kfolds_2.process_data()
kfolds_2.split(folds=5,seed=0)
for train, test in kfolds.return_all_splits():
    # Put your classifier here
    pass