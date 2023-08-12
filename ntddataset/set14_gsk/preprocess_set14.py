# -*- coding: utf-8 -*-
import os
import string
import pandas as pd
from rdkit import Chem
import numpy as np

datafile = 'KinetoBoxes_final.xlsx'
df = pd.read_excel(datafile,sheet_name=3,skiprows=1)

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

df = df.drop(list_to_remove, axis=0)

df['dataset'] = 'set14'
df['protein target'] = ''
df['treatment concentration'] = ''
df['treatment unit'] = ''
df['exposure time (h)'] = ''
df['type'] = ''
df['model type'] = ''
df['description'] = ''


chemical_info =      df[['SMILES',      'active ingredient','salt','protein target','treatment concentration','treatment unit','exposure time (h)','type','model type','description', 'dataset', 'ChEMBL_CMPD_NUMBER', 'BOX (2)', 'IFI', 'MW', 'MF', 'aring parent', 'clogp parent', 'hba parent', 'hbd parent', 'heavy parent', 'tpsa parent']]
chemical_info.columns = ['smiles code', 'active ingredient','salt','protein target','treatment concentration','treatment unit','exposure time (h)','type','model type','description' ,'dataset', 'ChEMBL_CMPD_NUMBER','BOX','IFI', 'MW', 'MF','aring parent', 'clogp parent', 'hba parent', 'hbd parent', 'heavy parent', 'tpsa parent']

#Anti - Trypanosoma cruzi activity

assay1 = chemical_info.copy(deep=True) 
assay1['assay value'] = df['pIC50  T. cruzi β-GAL']
assay1['measurement'] = 'pIC50  T. cruzi β-GAL'
assay1['target disease'] = 'Anti - Trypanosoma cruzi activity'
assay1['relation'] = df['modifier (1)']

assay2 = chemical_info.copy(deep=True) 
assay2['assay value'] = df['pIC50  T. cruzi imag INF CELL']
assay2['measurement'] = 'pIC50  T. cruzi imag INF CELL'
assay2['target disease'] = 'Anti - Trypanosoma cruzi activity'
assay2['relation'] = df['modifier (2)']

assay3 = chemical_info.copy(deep=True) 
assay3['assay value'] = df['pIC50  T. cruzi imag AM CELL']
assay3['measurement'] = 'pIC50  T. cruzi imag AM CELL'
assay3['target disease'] = 'Anti - Trypanosoma cruzi activity'
assay3['relation'] = df['modifier (3)']

assay4 = chemical_info.copy(deep=True) 
assay4['assay value'] = df['pIC50  T. cruzi imag CELL']
assay4['measurement'] = 'pIC50  T. cruzi imag CELL'
assay4['target disease'] = 'Anti - Trypanosoma cruzi activity'
assay4['relation'] = df['modifier (4)']

#Anti - Leishmania donovani activity

assay5 = chemical_info.copy(deep=True) 
assay5['assay value'] = df['pIC50 L. donovani FLINT']
assay5['measurement'] = 'pIC50 L. donovani FLINT'
assay5['target disease'] = 'Anti - Leishmania donovani activity'
assay5['relation'] = df['modifier (5)']

assay6 = chemical_info.copy(deep=True) 
assay6['assay value'] = df['pIC50 L. donovani Imag  INF MAC']
assay6['measurement'] = 'pIC50 L. donovani Imag  INF MAC'
assay6['target disease'] = 'Anti - Leishmania donovani activity'
assay6['relation'] = df['modifier (6)']

assay7 = chemical_info.copy(deep=True) 
assay7['assay value'] = df['pIC50 L. donovani Imag  AM MAC']
assay7['measurement'] = 'pIC50 L. donovani Imag  AM MAC'
assay7['target disease'] = 'Anti - Leishmania donovani activity'
assay7['relation'] = df['modifier (7)']

assay8 = chemical_info.copy(deep=True) 
assay8['assay value'] = df['pIC50 L. donovani Imag MAC']
assay8['measurement'] = 'pIC50 L. donovani Imag MAC'
assay8['target disease'] = 'Anti - Leishmania donovani activity'
assay8['relation'] = df['modifier (8)']

#anti - Trypanosoma brucei activity

assay9 = chemical_info.copy(deep=True) 
assay9['assay value'] = df['pIC50 T. brucei FLINT']
assay9['measurement'] = 'pIC50 T. brucei FLINT'
assay9['target disease'] = 'Anti - Trypanosoma brucei activity'
assay9['relation'] = df['modifier (9)']

assay10 = chemical_info.copy(deep=True) 
assay10['assay value'] = df['pIC50 T. brucei LUM']
assay10['measurement'] = 'pIC50 T. brucei LUM'
assay10['target disease'] = 'Anti - Trypanosoma brucei activity'
assay10['relation'] = df['modifier (10)']

#HepG2 toxicity

assay11 = chemical_info.copy(deep=True) 
assay11['assay value'] = df['pIC50 HepG2']
assay11['measurement'] = 'pIC50 HepG2'
assay11['target disease'] = 'HepG2 toxicity'
assay11['relation'] = df['modifier (11)']

#NIH 3T3 toxicity

assay12 = chemical_info.copy(deep=True) 
assay12['assay value'] = df['pIC50 3T3']
assay12['measurement'] = 'pIC50 3T3'
assay12['target disease'] = 'NIH 3T3 toxicity'
assay12['relation'] = df['modifier (13)']

#CYP51 activity

assay13 = chemical_info.copy(deep=True) 
assay13['assay value'] = df['pIC50 CYP51']
assay13['measurement'] = 'pIC50 CYP51'
assay13['target disease'] = 'CYP51 activity'
assay13['relation'] = df['modifier (12)']

df_final = pd.concat([assay1, assay2, assay3 , assay4 ,assay5 , assay6 ,assay7 , assay8 ,assay9 , assay10 ,assay11 , assay12 ,assay13], ignore_index=True)

chemical_info.columns = ['smiles code', 'active ingredient','salt','protein target','treatment concentration','treatment unit','exposure time (h)','type','model type','description' ,'dataset', 'ChEMBL_CMPD_NUMBER','BOX','IFI', 'MW', 'MF','aring parent', 'clogp parent', 'hba parent', 'hbd parent', 'heavy parent', 'tpsa parent']
df_final = df_final[['smiles code','active ingredient','salt','relation','assay value','target disease','protein target','treatment concentration','treatment unit','measurement','type','model type','description','dataset','ChEMBL_CMPD_NUMBER','BOX','IFI', 'MW', 'MF','aring parent', 'clogp parent', 'hba parent', 'hbd parent', 'heavy parent', 'tpsa parent']]

desc_df = pd.read_excel(datafile,sheet_name=4,skipfooter=4)

before = ''

list_to_remove = []

for i, row in desc_df.iterrows():        
    if row['COLUMN NAME'] == 'pIC50 T. cruzi CYP51':
        desc_df['COLUMN NAME'][i] = 'pIC50 CYP51'
    if pd.isnull(row['CHARACTERISTICS']):
        row['CHARACTERISTICS'] = before
    else:
        before = row['CHARACTERISTICS'] 

for i, row in df_final.iterrows():    
    try:        
        if row['relation'] != '=':
          list_to_remove.append(i)
          continue        
        desc = desc_df[((desc_df['CHARACTERISTICS'] ==  row['target disease']) & (desc_df['COLUMN NAME'] == row['measurement'])) ]
        description =  desc.iloc[0,2]
        method = desc.iloc[0,3]
        descStr = description
        if method != 'na':
            descStr = descStr + ', ' +  method
        df_final['description'][i] = descStr 
    except:
        print(i)
        break

df_final = df_final.drop(list_to_remove, axis=0)

df_final.to_excel('set14_result2.xlsx', # directory and file name to write
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

