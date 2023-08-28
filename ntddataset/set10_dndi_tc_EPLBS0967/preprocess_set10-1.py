# -*- coding: utf-8 -*-
import os
import sys
import string
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
import numpy as np


infofile = 'ASSAY_TABS.txt' # 실험 정보 
datafile = 'ACTIVITY_DETAIL_TABS.txt' #실험값 
sdffile = 'COMPOUND.sdf' # 분자 구조 

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
df = pd.read_csv(datafile,delimiter='\t')

# read assay infos
info = pd.read_csv(infofile, delimiter='\t')


list_to_remove.clear()

for index, row in df.iterrows():    
    compoundId = row['COMPOUND_ID']
    assayId = row['ASSAY_ID']
            
    # relation check
    if row['STATUS'].strip() != '=':
        #print(row['STATUS'])
        list_to_remove.append(index)
        continue        
    
    # compound id 로  smiles 값 조회 
    compound = sdf.loc[sdf['COMPOUND_ID'] == compoundId]
    
    if compound is None or compound.empty:
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
        df.loc[index, 'REFERENCE_ID']= assay.iloc[0,1]
        df.loc[index, 'ASSAY_DESCRIPTION']= assay.iloc[0,2]
        df.loc[index, 'ASSAY_TYPE']= assay.iloc[0,3] 
        df.loc[index, 'ASSAY_TEST_TYPE']= assay.iloc[0,4] 
        df.loc[index, 'ASSAY_ORGANISM']= assay.iloc[0,5]
        df.loc[index, 'ASSAY_STRAIN']= assay.iloc[0,6]
        df.loc[index, 'ASSAY_TAX_ID']= assay.iloc[0,7]
        df.loc[index, 'ASSAY_TISSUE']= assay.iloc[0,8]
        df.loc[index, 'ASSAY_CELL_LINE']= assay.iloc[0,9]
        df.loc[index, 'ASSAY_SUBCELLULAR_FRACTION']= assay.iloc[0,10]                
        df.loc[index, 'TARGET_TYPE']= assay.iloc[0,11]
        df.loc[index, 'TARGET_NAME']= assay.iloc[0,12]
        df.loc[index, 'TARGET_ACCESSION']= assay.iloc[0,13]
        df.loc[index, 'TARGET_ORGANISM']= assay.iloc[0,14]
        df.loc[index, 'TARGET_TAX_ID']= assay.iloc[0,15]
        df.loc[index, 'SRC_ID']= assay.iloc[0,16]
        df.loc[index, 'ASSAY_SOURCE']= assay.iloc[0,17]
           

df = df.drop(list_to_remove, axis=0)   

df['dataset'] = 'set10-1'
df['protein target'] = ''
df['treatment concentration'] = ''
df['exposure time (h)'] = ''
df['model type'] = ''

df = df[[     'smiles',      'active ingredient','salt','STATUS', 'VALUE',      'TARGET_NAME',    'protein target','treatment concentration','UNIT','exposure time (h)',  'MEASUREMENT_TYPE', 'ASSAY_TEST_TYPE','model type','ASSAY_DESCRIPTION', 'dataset','REFERENCE_ID','ASSAY_TYPE','ASSAY_ORGANISM','ASSAY_STRAIN','TARGET_TYPE','SRC_ID','ASSAY_SOURCE','UPPER_VALUE','SD_MINUS','SD_PLUS','ACTIVITY_COMMENT']]
df.columns = ['smiles code', 'active ingredient','salt','relation','assay value','target disease', 'protein target','treatment concentration','treatment unit','exposure time (h)','measurement','type','model type','description' ,'dataset','REFERENCE_ID','ASSAY_TYPE','ASSAY_ORGANISM','ASSAY_STRAIN','TARGET_TYPE','SRC_ID','ASSAY_SOURCE','UPPER_VALUE','SD_MINUS','SD_PLUS','ACTIVITY_COMMENT']

df.to_csv('set10-1_result.csv', index=False)

