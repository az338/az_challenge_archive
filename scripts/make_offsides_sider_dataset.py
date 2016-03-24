
import pandas as pd
DATA_DIR = '/scratch/az338/ucc-fileserver/side_effects_data/'

# load SE datasets

offsides = pd.read_table(DATA_DIR+'3003377s-offsides.tsv')
sider_se = pd.read_table(DATA_DIR+'meddra_all_se.tsv',header=None)
sider_se = sider_se.drop(0,1)
sider_se = sider_se[sider_se[3] == 'PT']

#sider_fq = pd.read_table('/home/az338/ucc-fileserver/sider/meddra_freq.tsv',header=None)
#sider_fq = sider_fq.drop(0,1)


# concatenate sider and offsides side effect
sider_se = sider_se.drop([3,4],1)
sider_se.columns = ['stitch_id','umls_id','event']
allSE = pd.concat([sider_se,offsides[['stitch_id','umls_id','event']]])
allSE.to_csv(DATA_DIR+'sider_and_offsides_SE.csv')

# load stitch_id  smiles mapping dataset 
stitch_smiles = pd.read_table('/home/az338/ucc-fileserver/stitch_v4/chemicals.v4.0.tsv')


# match stitch id in SE dataset to smiles
idx_offsides = pd.match(offsides['stitch_id'],stitch_smiles['chemical'])
idx_offsides = idx_offsides[idx_offsides >= 0]

idx_sider = pd.match(sider_se['stitch_id'],stitch_smiles['chemical'])
idx_sider = idx_sider[idx_sider >= 0]

# concatenate sider and offsides structures
offsides_struct = stitch_smiles.iloc[idx_offsides][['chemical','SMILES_string']].drop_duplicates()
sider_struct = stitch_smiles.iloc[idx_sider][['chemical','SMILES_string']].drop_duplicates()
all_struct = pd.concat([offsides_struct,sider_struct]).drop_duplicates()
all_struct.to_csv(DATA_DIR+'sider_offsides_smiles.csv')


