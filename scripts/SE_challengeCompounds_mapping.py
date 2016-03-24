import pandas as pd

DATA_DIR = '/home/az338/ucc-fileserver/AZ_challenge_data/'
SE_DATA_DIR = '/scratch/az338/ucc-fileserver/side_effects_data/'

# load target/SE mapping 
target_SE_map= pd.read_csv(DATA_DIR+'target_side_effects.csv')

# load challenge compound-target mapping

challenge_cmp_target_map = pd.read_csv(DATA_DIR+'drug_synergy_data/Drug_info_release_curated.csv')
challenge_cmp_target_map.columns = ['ChallengeName','Target'] + list(challenge_cmp_target_map.columns[2:])
# make mapping flat (separate targets to different lines)
cmp_target_map = challenge_cmp_target_map[['ChallengeName','Target']]

flat_map = []
for i, r in cmp_target_map.iterrows():
    for t in r['Target'].split(','):
        flat_map.append([r['ChallengeName'],t.rstrip(' ').lstrip(' ')])

flat_map = pd.DataFrame(flat_map).drop_duplicates()
flat_map.columns = ['ChallengeName','Target']

# takes a data frame of SE profiles for the targets
# of a given compound and merge them (using OR operation) ` 
def merge_target_SEprofiles(df) :
    idx = pd.match(df['Target'],target_SE_map['symbol'])
    idx = idx[idx >= 0]
    if len(idx) == 0 :
        return 
    SE = target_SE_map.iloc[idx]
    SE = SE.drop('symbol',axis=1)
    SE = SE.values.sum(axis = 0) # sum in columns
    SE[SE > 0] = 1  # consider everything > 0 as 1 (OR operation)
    return [df['ChallengeName'].values[0]]+SE.tolist()

cmp_SE_map = pd.DataFrame(flat_map.groupby('ChallengeName').apply(merge_target_SEprofiles).dropna().tolist())
cmp_SE_map.columns = target_SE_map.columns

cmp_SE_map.to_csv(DATA_DIR+'challenge_compounds_SE.csv')
