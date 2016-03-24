
import pandas as pd
import MySQLdb
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


DATA_DIR = '/home/az338/ucc-fileserver/AZ_challenge_data/'
SE_DATA_DIR = '/scratch/az338/ucc-fileserver/side_effects_data/'

# load compounds from side effect dataset
SE_data = pd.read_csv(SE_DATA_DIR+'/sider_and_offsides_SE.csv')
SE_struct = pd.read_csv(SE_DATA_DIR+'sider_offsides_smiles.csv')

# compute ECFP6 fingerprints for molecules of the SE dataset
se_mols = [Chem.MolFromSmiles(s) for s in SE_struct['SMILES_string']]
se_fp = [AllChem.GetMorganFingerprint(m,3) for m in se_mols]
se_id = SE_struct['chemical']

# load uniprot for overlapping challenge targets in chembl20
targets = pd.read_table(DATA_DIR+'chembl20_challenge_overlaping.txt',sep=' ')


# connect to chembl20
conn = MySQLdb.connect("localhost","az338","pleasechange","chembl_20")
c = conn.cursor()

compound_target_mapping = []
for unip,tname in zip(targets['UNIPROT'].tolist(),targets['SYMBOL']):
    # map uniprots to tids
    c.execute('select td.tid, td.pref_name, accession from component_sequences as cs, target_components as tc, target_dictionary as td where cs.component_id=tc.component_id and tc.tid=td.tid and accession="'+unip+'"')
    tid_map = map(list,c.fetchall())

    # map tid to assay ids
    aid_map = []
    for tid in [l[0] for l in tid_map] :
        c.execute('select assay_id, tid, description, assay_type  from assays where tid="'+str(tid)+'" and assay_type in ("B","F")')
        aid_map.append(map(list,c.fetchall()))
    aid_map = sum(aid_map,[]) # again unnest the results

    # map assay ids to active mol reg no
    mol_map = []
    for aid in [l[0] for l in aid_map] :
        c.execute('select molregno, assay_id, concat(standard_relation, round(standard_value,2)), standard_type, standard_units from activities where pchembl_value>5 and standard_type in ("EC50","IC50","Ki","Kd") and assay_id="'+str(aid)+'"')
        mol_map.append(map(list,c.fetchall()))
    mol_map = sum(mol_map,[])
    
    # get unique mol reg no.
    molregno = list(set([l[0] for l in mol_map]))

    # get SMILES/structure corresponding to those mol reg no
    struct = []
    for m in molregno:
        c.execute('select molregno, canonical_smiles from compound_structures where molregno="'+str(m)+'"')
        struct.append(map(list,c.fetchall())) 
    struct = sum(struct,[])

    # compute fingerprints
    mols = [Chem.MolFromSmiles(s[1]) for s in struct]
    fp = [AllChem.GetMorganFingerprint(m,3) for m in mols]

    # for each compound in side effect dataset compute similarity
    # to target known actives 
    sim = [DataStructs.BulkTanimotoSimilarity(ch,fp) for ch in se_fp]


    # store number of target known active for each compound in side effect dataset
    for i in xrange(len(sim)):
        cn = 0
        for s in sim[i]:
            if s >= .6 : cn+=1
        if cn > 0:
            compound_target_mapping.append([se_id[i],unip,tname,cn])

# convert compound/target mapping to data frame
compound_target_mapping = pd.DataFrame(compound_target_mapping)
compound_target_mapping.columns = ['stitch_id','uniprot','symbol','n']

# filter data to only compound similar to more than 5 
filtered_mapping = compound_target_mapping[compound_target_mapping['n'] >= 1]

# target to SE mapping 
target_se_map = pd.merge(filtered_mapping,SE_data,on='stitch_id',how='inner')
target_se_map['val'] = 1
target_se_mat = pd.pivot_table(target_se_map, rows='symbol', values='val', cols='event', fill_value=0)
target_se_mat.to_csv(DATA_DIR+'target_side_effects.csv')





