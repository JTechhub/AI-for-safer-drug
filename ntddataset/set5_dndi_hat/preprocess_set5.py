# -*- coding: utf-8 -*-
import os
import sys
import string
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
import numpy as np


infofile = 'Assays_table.xlsx'
datefile = 'Compounds_table.txt'
sdffile = 'Structures_Table.sdf'

# read sdf
sdf = PandasTools.LoadSDF(sdffile)
list_to_remove = []

for index, row in sdf.iterrows():    
    mol = row['ROMol']
    smi = Chem.MolToSmiles(mol)
    sdf.loc[index, 'smiles']=smi
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

            sdf.loc[index, 'active ingredient']=active
            sdf.loc[index, 'salt']=salt

        except:
            print (index, smi)
            list_to_remove.append(index)
    else:
        sdf.loc[index, 'active ingredient']=smi
        sdf.loc[index, 'salt']='.'
        
sdf = sdf.drop(list_to_remove, axis=0)


# read data file
df = pd.read_csv(datefile,sep='\t')

# read assay infos
info = pd.read_excel(infofile,sheet_name=0)

list_to_remove.clear()

for index, row in df.iterrows():    
    compoundId = row['COMPOUND_ID']
    assayId = row['ASSAY_ID']
            
    # relation check
    if row['STATUS'] != '=':
        list_to_remove.append(index)
        continue        
    
    # compound id 로  smiles 값 조회 
    compound = sdf.loc[sdf['COMPOUND_ID'] == compoundId]
    
    if compound is None:
        print ("no compound",index, compoundId)
        list_to_remove.append(index)
        continue
    else:    
        df.loc[index, 'active ingredient']= compound.iloc[0,4]
        df.loc[index, 'salt']= compound.iloc[0,5]
        df.loc[index, 'smiles']= compound.iloc[0,3]
        
    # assay id 로  실험 정보 조회 
    assay = info.loc[info['ASSAY_ID'] == assayId]
    
    if assay is None:
        print ("no assay",index,assayId)
        list_to_remove.append(index)
        continue        
    else:    
        df.loc[index, 'ASSAY_DESCRIPTION']= assay.iloc[0,1]
        df.loc[index, 'ASSAY_TYPE']= assay.iloc[0,2] 
        df.loc[index, 'ASSAY_ORGANISM']= assay.iloc[0,3]
        df.loc[index, 'ASSAY_STRAIN']= assay.iloc[0,4]
        df.loc[index, 'TARGET_TYPE']= assay.iloc[0,5]
        df.loc[index, 'TARGET_NAME']= assay.iloc[0,6]
        df.loc[index, 'MEASUREMENT_TYPE']= assay.iloc[0,7]
        df.loc[index, 'DOSE']= assay.iloc[0,8]
        df.loc[index, 'dose_unit']= assay.iloc[0,9]
        df.loc[index, 'TIMEPOINT']= assay.iloc[0,10]
        df.loc[index, 'timepoint_unit']= assay.iloc[0,11]
        df.loc[index, 'ASSAY_SOP']= assay.iloc[0,13]
        df.loc[index, 'New ASSAY_SOP']= assay.iloc[0,14]
        df.loc[index, 'Comment']= assay.iloc[0,15]
    
            
df = df.drop(list_to_remove, axis=0)   

df['dataset'] = 'set5'
df['protein target'] = ''
df['treatment concentration'] = ''
df['exposure time (h)'] = ''
df['type'] = ''
df['model type'] = ''

df = df[[     'smiles',      'active ingredient','salt','STATUS', 'VALUE',      'TARGET_NAME',    'protein target','treatment concentration','UNIT','exposure time (h)',  'MEASUREMENT_TYPE', 'type','model type','ASSAY_DESCRIPTION', 'dataset','ASSAY_TYPE','ASSAY_STRAIN','TARGET_TYPE','DOSE','timepoint_unit','SUPPLIER','ASSAY_SOP','New ASSAY_SOP','Comment']]
df.columns = ['smiles code', 'active ingredient','salt','relation','assay value','target disease', 'protein target','treatment concentration','treatment unit','exposure time (h)','measurement','type','model type','description' ,'dataset','ASSAY_TYPE','ASSAY_STRAIN','TARGET_TYPE','DOSE','timepoint_unit','SUPPLIER','ASSAY_SOP','New ASSAY_SOP','Comment']


df.to_excel('dataset5.xlsx', # directory and file name to write
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


