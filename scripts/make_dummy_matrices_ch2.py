
# create dummy leaderboard and prediciton matrix for challenge2
import pandas as pd
import numpy as np
import numpy.random as npr
from copy import copy

# load leaderboard file (observed synergies)
lb_ch2 = pd.read_csv('/home/az338/Desktop/Link to ucc-fileserver/AZ_challenge_data/drug_synergy_data/ch2_leaderBoard_monoTherapy.csv')
cell_lines = lb_ch2['CELL_LINE'].drop_duplicates().values
combinations = lb_ch2['COMBINATION_ID'].drop_duplicates().values

# make random synergy prediction matrix
rand_mat = npr.rand(len(combinations),len(cell_lines))
syn_mat = copy(rand_mat)
# convert to binary prediction
syn_mat[syn_mat < .5] = 0
syn_mat[syn_mat >= .5] = 1

# put columns and row names into synergy matrix
syn_mat = pd.DataFrame(syn_mat)
syn_mat.columns = cell_lines
syn_mat.index = combinations

# put columns and row names into 
conf_mat = pd.DataFrame(rand_mat)
conf_mat.columns = cell_lines
conf_mat.index = combinations

# generate random scores for the synergies (between 80 and -50)
rand_scores = npr.rand(len(lb_ch2))*(80+50)-50
lb_ch2['SYNERGY_SCORE'] = rand_scores

# export to csv
lb_ch2.to_csv('/home/az338/ucc-fileserver/AZ_challenge_data/drug_synergy_data/ch2_leaderBoard_dummy.csv')
syn_mat.to_csv('/home/az338/ucc-fileserver/AZ_challenge_data/ch2_pred_dummy.csv')
conf_mat.to_csv('/home/az338/ucc-fileserver/AZ_challenge_data/ch2_confidence_dummy.csv')
