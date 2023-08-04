import os
import pandas as pd

from rdkit import Chem

wdir = r'directory/for/ChEMBLNTD/set1_gsk'
csv = 'chemblntd_all.txt'
df = pd.read_csv(os.path.join(wdir, csv), delimiter='\t')

#Extracting smiles code, active ingredient, salt
list_to_remove = []
for index, row in df.iterrows():
    smi = row['SMILES']
    mol = Chem.MolFromSmiles(smi)
    if '.' in smi:
        try:
            fragments_list = list(Chem.GetMolFrags(mol, asMols=True))
            fragments_sizes = [frag.GetNumAtoms() for frag in fragments_list]

            largest_i = fragments_sizes.index(max(fragments_sizes))
            largest_frag = fragments_list[largest_i]

            active = Chem.MolToSmiles(largest_frag)

            fragments_list.remove(largest_frag)
            salt_smi_list = [Chem.MolToSmiles(frag) for frag in fragments_list]
            salt = '.'.join(salt_smi_list)
            #print(index, smi, active, salt)

            df.loc[index, 'active ingredient']=active
            df.loc[index, 'salt']=salt
        except:
            print (index, smi)
            list_to_remove.append(index)
    else:
        df.loc[index, 'active ingredient']=smi
        df.loc[index, 'salt']='.'

#remove smiles codes that cause error in rdkit.
df = df.drop(list_to_remove, axis=0)

df.to_csv(os.path.join('/directory/for/ChEMBLNTD_curation','set1_raw.csv'), index=False)
