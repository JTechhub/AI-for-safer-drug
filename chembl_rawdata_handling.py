# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:14:56 2022

@author: Gil
"""
import os
import pandas as pd
from rdkit import Chem

#wdir is short for workding directory.
wdir = r'/directory/where/csv/file/is/'
csv = 'DOWNLOAD-n0i9yh99WaQDBaNwN8GG-_4iqK6LwkZOMsWMMKaj2BI=.csv'

"""
Error occured since pd.read_csv assumes that csv file is separated by ','.
If delimiter is not ',', it should be specifically inputted to 'sep' in pd.read_csv
"""
#df = pd.read_csv(os.path.join(wdir, csv))

df = pd.read_csv(os.path.join(wdir, csv), sep=';')
print ('Initial data size:', df.shape)
print ('*'*40)

#To check column names of csv files
print ('Check the list of column names in the csv file:', df.columns)
print ('*'*40)

#Data Validity Comment is empty if data is normal. 
print ('Statistics of Data Validity Comment\n:',df['Data Validity Comment'].value_counts())
print ('*'*40)

#Standard Relation contains =, <, >, and so forth.
print ('Statistics of Standard Relation\n',df['Standard Relation'].value_counts())
print ('*'*40)

"""
The column the value is '=', not just =.
Therefore matching pattern should be '=', not =
"""
#df = df[df['Standard Relation']=='=']
df = df[df['Standard Relation']=='\'=\'']
df = df[df['Data Validity Comment'].isna()]

print ('Statistics of unit\n',df['Standard Units'].value_counts())
print ('*'*40)

#Two different units were fouond: ug/mL and nM.
#Unit should be standardized in order to compare 'Standard Value'.
df_diffunit = df[df['Standard Units']=='ug.mL-1']
print (df_diffunit.shape)

df_nM = df[df['Standard Units']!='ug.mL-1']
print (df_nM.shape)
print ('*'*40)

"""
Calculation failed with the commented codes below since the type was str.
"""
#converted = 10**(6)*df_diffunit['Standard Value']/df_diffunit['Molecular Weight']
#print (df_diffunit['Standard Value'].to_list(), type(df_diffunit['Standard Value'].to_list()[0]))
#print (df_diffunit['Molecular Weight'].to_list(), type(df_diffunit['Molecular Weight'].to_list()[0]))

converted = 10**(6)*df_diffunit['Standard Value']/df_diffunit['Molecular Weight'].astype('float') #Equation for unit conversion

print ('Maximum IC50 value:',df['Standard Value'].max())
print ('Minimum IC50 value:',df['Standard Value'].min())

df_diffunit['Standard Value']=converted
df_diffunit['Standard Units']='nM'
print (df_diffunit[['Standard Value','Standard Units']])
print ('*'*40)
df = pd.concat([df_nM, df_diffunit])

#In this work, our goal is to identify active ingredient for the drug.
#That activie ingredient should be a single molecule.
#Some of smiles code may contain more than one molecular structure.
#In that case, smiels code contains period, which is a sign that the molecules are not connected.
for index, row in df[df['Smiles'].str.contains('\.')].iterrows():
    smi_list = row['Smiles'].split('.')
    
    active_ingredient = smi_list[0]
    for smi in smi_list[1:]:
        if len(active_ingredient) < len(smi):
            active_ingredient = smi
    df.loc[index,'Smiles'] = active_ingredient

#Check if metal complexes in the data.
#Normally, druggable chemical space only includes atoms below.
organic_set = ['H','C','N','O','F','P','S','Cl','Br','I']

df['metal']=0
for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['Smiles'])
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in organic_set:
            df.loc[index,'metal']=1
            break

metal_flag = df['metal'].sum()
if metal_flag == 0:
    print ('The number of molecules containing metal atom:', df['metal'].sum())
else:
    print ('The number of metal complexes:', df.loc[df['metal']!=0, 'Smiles'])
    df = df[df['metal'] == 0]
    print ('These are removed.')

#Let's check duplicated structures.
print ('ChEMBL ID duplicated data:', df.duplicated(subset='Molecule ChEMBL ID',keep=False).sum())
print ('Smiles duplicated data:', df.duplicated(subset='Smiles', keep=False).sum())

#Simple way to identify same structures.
"""
#One molecular structure can be translated into more than one smiles codes.
smi_list = ['C1=NC2=C(N1)C(=S)N=CN2',
            'n1c[nH]c2c1NC=NC2=S',
            '[H]N1C=NC2=C1C(=S)N=C([H])N2[H]']
for smi in smi_list:
    m = Chem.MolFromSmiles(smi)
    print (Chem.MolToSmiles(m))
"""

#This code regenerate smiles code using rdkit for smiles code standardization.
#RDKit has its own way to convert molecular structure into smiles code.
#By generating smiles code with rdkit Chem.MolToSmiles, standardized smiles code is generated.
for index, row in df.iterrows():
    m = Chem.MolFromSmiles(row['Smiles'])
    df.loc[index, 'rdkit_smi']=Chem.MolToSmiles(m)
print ('standardized smiles code duplicated data:', df.duplicated(subset='rdkit_smi',keep=False).sum())

#Check if ChEMBL ID, smiles code, and rdkit standardized smiels code can find identical number of duplicated data.
print ('Duplicated compounds between Smiles and rdkit_smi',(df[df.duplicated(subset='rdkit_smi',keep=False)].index==df[df.duplicated(subset='Smiles',keep=False)].index).sum())

import numpy as np
#pIC50 converstion
df['pIC50'] = -np.log10(df['Standard Value'])

#Check duplicated structure. 
df_dup = df[df.duplicated(subset='Smiles',keep=False)]
df_single = df[~df.duplicated(subset='Smiles',keep=False)]

print ('Check an example of duplicated structures.')
dup_smi_list = list(set(df_dup['Smiles'].to_list()))
print (df_dup[df_dup['Smiles']==dup_smi_list[0]])
print (df_dup.loc[df_dup['Smiles']==dup_smi_list[0],'Standard Value'])
print (df_dup.loc[df_dup['Smiles']==dup_smi_list[0],'pIC50'])

df_list = []
for smi in dup_smi_list:
    #Take dataframe with identical smiles code.
    temp_df = df_dup[df_dup['Smiles']==smi]
    
    #Calculate mean, max, min value of pIC50
    mean_pIC50 = -np.log10(temp_df['Standard Value'].mean())
    max_pIC50 = temp_df['pIC50'].max()
    min_pIC50 = temp_df['pIC50'].min()
    
    dif_mean_max = abs(max_pIC50 - mean_pIC50)
    dif_mean_min = abs(mean_pIC50 - min_pIC50)
    if max(dif_mean_max, dif_mean_min) < 0.5:
        #Get the IC50 mean value.
        temp_df['pIC50']=mean_pIC50
        #Collect dataframe.
        df_list.append(temp_df.iloc[0])

#Concatenate dataframes in the list.        
df_dup_final = pd.concat(df_list, axis=1).T

df_final = pd.concat([df_dup_final, df_single])

#Final data
print ('Final data size:', df_final.shape)
df_final.to_csv(os.path.join(wdir, 'preprocessed_'+csv),index=False)