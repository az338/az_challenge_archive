import pandas as pd

DATA_DIR = '/scratch/az338/ucc-fileserver/AZ_challenge_data/'

# challenge training data, cell/disease area(DA) mapping and challenge cmp/target mapping
challenge_pairs = pd.read_csv(DATA_DIR+'drug_synergy_data/ch1_train_combination_and_monoTherapy.csv')
cell_da_map = pd.read_csv(DATA_DIR+'sanger_molecular_data/cell_info.csv')
challenge_cmp_target_map = pd.read_csv(DATA_DIR+'drug_synergy_data/Drug_info_release_curated.csv')
challenge_cmp_target_map.columns = ['ChallengeName','Target'] + list(challenge_cmp_target_map.columns[2:])

# map disease area to corresponding cell-line
challenge_pairs['DISEASE_AREA'] = cell_da_map.iloc[pd.match(challenge_pairs['CELL_LINE'], cell_da_map['Sanger.Name'])]['Disease.Area'].tolist()

# make mapping flat (separate targets to different lines)
cmp_target_map = challenge_cmp_target_map[['ChallengeName','Target']]

flat_map = []
for i, r in cmp_target_map.iterrows():
    for t in r['Target'].split(','):
        flat_map.append([r['ChallengeName'],t.rstrip(' ').lstrip(' ')])

flat_map = pd.DataFrame(flat_map).drop_duplicates()
flat_map.columns = ['Compound','Target']

# convert compound-compound associations to target-target associations
synergy_scores = challenge_pairs[['DISEASE_AREA','COMPOUND_A','COMPOUND_B','SYNERGY_SCORE','QA','CELL_LINE']]
target_synergy = []
for i, r in synergy_scores.iterrows():
    #targets_A = flat_map.iloc[pd.match(r['COMPOUND_A'],flat_map['Compound'])]['Target']
    #targets_B = flat_map.iloc[pd.match(r['COMPOUND_B'],flat_map['Compound'])]['Target']
    targets_A = flat_map[flat_map['Compound'] == r['COMPOUND_A']]['Target']
    targets_B = flat_map[flat_map['Compound'] == r['COMPOUND_B']]['Target']
    for tA in targets_A :
        for tB in targets_B :
            if r['QA']== 1:
                target_synergy.append([r['DISEASE_AREA'],tA,tB,r['SYNERGY_SCORE'],r['CELL_LINE']])

target_synergy = pd.DataFrame(target_synergy)
target_synergy.columns = ['DISEASE_AREA','TARGET_A','TARGET_B','SYNERGY','CELL_LINE']
target_synergy.to_csv(DATA_DIR+'target_synergy_train_by_disease_area.csv',index=False)

# average synergy scores per DA and target pair
#target_synergy = target_synergy.groupby(['DISEASE_AREA','TARGET_A','TARGET_B']).mean().reset_index()

# apply cutoff on which to consider a target pair synergistic or not
# uncomment the two lines below if binarization required
#target_synergy.ix[target_synergy['SYNERGY'] < 10,'SYNERGY'] = 0
#target_synergy.ix[target_synergy['SYNERGY'] >= 10,'SYNERGY'] = 1 

# expand synergy to synergy.disease area variables
#pivot =  target_synergy.pivot_table(rows=['TARGET_A','TARGET_B'], cols='DISEASE_AREA', values='SYNERGY').reset_index()
#pivot.columns = ['Target_A','Target_B'] + list(pivot.columns[2:]+'.Synergy')
#pivot.to_csv(DATA_DIR+'target_synergy_train_by_disease_area.csv',index=False)



