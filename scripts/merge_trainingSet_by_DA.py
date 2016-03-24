import pandas as pd

DATA_DIR = '/scratch/az338/ucc-fileserver/AZ_challenge_data/'

# challenge training data, cell/disease area(DA) mapping and challenge cmp/target mapping
challenge_pairs = pd.read_csv(DATA_DIR+'drug_synergy_data/ch1_train_combination_and_monoTherapy.csv')
cell_da_map = pd.read_csv(DATA_DIR+'sanger_molecular_data/cell_info.csv')

# map disease area to corresponding cell-line
challenge_pairs['DISEASE_AREA'] = cell_da_map.iloc[pd.match(challenge_pairs['CELL_LINE'], cell_da_map['Sanger.Name'])]['Disease.Area'].tolist()

# convert compound-compound associations to target-target associations
synergy_scores = challenge_pairs[['DISEASE_AREA','COMPOUND_A','COMPOUND_B','SYNERGY_SCORE','QA','CELL_LINE']]
# remove combinations with low quality 
synergy_scores  = synergy_scores[synergy_scores['QA'] == 1]


# write to file
synergy_scores.to_csv(DATA_DIR+'compound_synergy_trainSet_by_disease_area.csv',index=False)

# apply cutoff on which to consider a target pair synergistic or not
# uncomment the two lines below if binarization required
#synergy_scores.ix[target_synergy['SYNERGY'] < 10,'SYNERGY'] = 0
#synergy_scores.ix[target_synergy['SYNERGY'] >= 10,'SYNERGY'] = 1 

# expand synergy to synergy.disease area variables (TO REMOVE???)
#pivot =  synergy_scores.pivot_table(rows=['COMPOUND_A','COMPOUND_B'], cols='DISEASE_AREA', values='SYNERGY_SCORE').reset_index()
#pivot.columns = ['Compound_A','Compound_B'] + list(pivot.columns[2:]+'.Synergy')
#pivot.to_csv(DATA_DIR+'target_synergy_train_by_disease_area.csv',index=False)



