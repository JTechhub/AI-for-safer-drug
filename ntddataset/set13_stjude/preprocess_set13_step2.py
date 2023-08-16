# -*- coding: utf-8 -*-
import os
import string
import pandas as pd
from rdkit import Chem
import numpy as np

datafile = 'dataset13.xlsx'
df = pd.read_excel(datafile,sheet_name=0)

list_to_remove = []

for index, row in df.iterrows():    
    if row['Relation'] != '=':
      list_to_remove.append(index)
      continue            
    
    smi = row['MOLSMILES']
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

df = df.drop(list_to_remove, axis=0)


df['dataset'] = 'set13'
df['protein target'] = ''
df['treatment concentration'] = ''
df['exposure time (h)'] = ''
df['type'] = ''
df['model type'] = ''


df = df[['MOLSMILES',  'active ingredient','salt','Relation','Standard Value','Target Organism', 'protein target','treatment concentration','Standard Units','exposure time (h)','Standard Type', 'type','model type','Description', 'dataset','Assay Type','Assay Src Description','Assay Organism','Target Type','Target Name','Reference','CMPD_ID','STJUDE_ID','conc_uM','Growth Inhibition(%)','DataSet','ID']]
df.columns = ['smiles code', 'active ingredient','salt','relation','assay value','target disease', 'protein target','treatment concentration','treatment unit','exposure time (h)','measurement','type','model type','description' ,'dataset','Assay Type','Assay Src Description','Assay Organism','Target Type','Target Name','Reference','CMPD_ID','STJUDE_ID','conc_uM','Growth Inhibition(%)','DataSet','ID']


df.to_excel('set13_result.xlsx', # directory and file name to write
            sheet_name = 'Sheet1', 
            na_rep = '', 
            float_format = "%.2f", 
            header = True, 
            #columns = ["group", "value_1", "value_2"], # if header is False
            index = False, 
            #index_label = "id", 
            startrow = 0, 
            startcol = 0, 
            #engine = 'xlsxwriter', 
            #freeze_panes = (2, 0)
            ) 

