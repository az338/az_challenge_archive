import pandas as pd
import numpy as np

DATA_DIR = '/home/az338/ucc-fileserver/AZ_challenge_data/'

# load side effect profiles of individual drugs
se_matrix = pd.read_csv(DATA_DIR+'challenge_compounds_SE.csv')

# load all challenge training and test sets
train_ch1 = pd.read_csv(DATA_DIR+'drug_synergy_data/ch1_train_combination_and_monoTherapy.csv')
train_ch1 = train_ch1[train_ch1['QA'] == 1]
leaderboard_ch1 = pd.read_csv(DATA_DIR+'drug_synergy_data/ch1_leaderBoard_monoTherapy_dummy.csv')
leaderboard_ch2 = pd.read_csv(DATA_DIR+'drug_synergy_data/ch2_leaderBoard_monoTherapy.csv')

# get unique combination ids
combination_ids_tr1 = np.unique([train_ch1.COMBINATION_ID.values+'>'+train_ch1.CELL_LINE.values])
combination_ids_l1 =  np.unique([leaderboard_ch1.COMBINATION_ID.values+'>'+leaderboard_ch1.CELL_LINE.values])
combination_ids_l2 =  np.unique([leaderboard_ch2.COMBINATION_ID.values+'>'+leaderboard_ch2.CELL_LINE.values])

# get number of SE/features (n: 1022)
n_features = len(se_matrix.columns)-1


# combine side effect profiles for each compound in the combination id list 
combi_SE_tr1 = []
for c in combination_ids_tr1 :
    c1 = c.replace('>','.').split('.')[0] # first compound of the combination
    c2 = c.replace('>','.').split('.')[1]  # second compound of the combination
    SE1idx  = pd.match([c1],se_matrix['symbol'].values)
    SE2idx = pd.match([c2],se_matrix['symbol'].values)
    if SE1idx >= 0 :
        SE1 = se_matrix.iloc[SE1idx].drop('symbol',axis=1).values
    else :
        SE1 = np.array([0]*n_features)
    if SE2idx >= 0 :
        SE2 = se_matrix.iloc[SE2idx].drop('symbol',axis=1).values
    else :
        SE2 = np.array([0]*n_features)
    combi_SE_tr1.append(np.append(c,SE1 + SE2))

combi_SE_l1 = []
for c in combination_ids_l1 :
    c1 = c.replace('>','.').split('.')[0] # first compound of the combination
    c2 = c.replace('>','.').split('.')[1]  # second compound of the combination
    SE1idx  = pd.match([c1],se_matrix['symbol'].values)
    SE2idx = pd.match([c2],se_matrix['symbol'].values)
    if SE1idx >= 0 :
        SE1 = se_matrix.iloc[SE1idx].drop('symbol',axis=1).values
    else :
        SE1 = np.array([0]*n_features)
    if SE2idx >= 0 :
        SE2 = se_matrix.iloc[SE2idx].drop('symbol',axis=1).values
    else :
        SE2 = np.array([0]*n_features)
    combi_SE_l1.append(np.append(c,SE1 + SE2))

combi_SE_l2 = []
for c in combination_ids_l2 :
    c1 = c.replace('>','.').split('.')[0] # first compound of the combination
    c2 = c.replace('>','.').split('.')[1]  # second compound of the combination
    SE1idx  = pd.match([c1],se_matrix['symbol'].values)
    SE2idx = pd.match([c2],se_matrix['symbol'].values)
    if SE1idx >= 0 :
        SE1 = se_matrix.iloc[SE1idx].drop('symbol',axis=1).values
    else :
        SE1 = np.array([0]*n_features)
    if SE2idx >= 0 :
        SE2 = se_matrix.iloc[SE2idx].drop('symbol',axis=1).values
    else :
        SE2 = np.array([0]*n_features)
    combi_SE_l2.append(np.append(c,SE1 + SE2))

combi_df_tr = pd.DataFrame(combi_SE_tr1)
combi_df_tr.columns = ['ID'] + se_matrix.columns[1:].tolist()
#combi_df.rename(columns={0:'COMBINATION_ID'},inplace=True)
combi_df_tr.to_csv(DATA_DIR + 'az_combination_train_side_effects.csv', index=False)

combi_df_t1 = pd.DataFrame(combi_SE_l1)
combi_df_t1.columns = ['ID'] + se_matrix.columns[1:].tolist()
combi_df_t1.to_csv(DATA_DIR + 'az_combination_test_1_side_effects.csv', index=False)

combi_df_t2 = pd.DataFrame(combi_SE_l2)
combi_df_t2.columns = ['ID'] + se_matrix.columns[1:].tolist()
combi_df_t2.to_csv(DATA_DIR + 'az_combination_test_2_side_effects.csv', index=False)
