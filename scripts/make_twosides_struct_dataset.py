import pandas as pd


DATA_DIR = '/home/az338/ucc-fileserver/AZ_challenge_data/'
SE_DATA_DIR = '/scratch/az338/ucc-fileserver/side_effects_data/'

# load twosides
SE_data = pd.read_table(SE_DATA_DIR+'3003377s-twosides.tsv')

# load compound structures 
SE_struct = pd.read_table('/home/az338/ucc-fileserver/stitch_v4/chemicals.v4.0.tsv')

# get unique compound identifiers in twosides
twosides_compounds = list(set(SE_data['stitch_id1'].values.tolist()+SE_data['stitch_id2'].values.tolist()))

# match these compounds to the OFFSIDES/SIDER structural dataset 
# to match the stitch IDs to their structure
idx = pd.match(twosides_compounds,SE_struct['chemical'])
twosides_struct = SE_struct.iloc[idx[idx>0]]
twosides_struct=twosides_struct.drop('molecular_weight',axis=1)

twosides_struct.to_csv(SE_DATA_DIR+'twosides_smiles.csv')
