# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 20:26:21 2019

@author: Gil
"""
import os

from rdkit_descriptor_calculation import *

wdir = r'/directory/where/csv/file/is/'
target = ['train.csv', 'test.csv']

for t in target:
    df = pd.read_csv(os.path.join(wdir,t))
    df_id = df[['Molecule ChEMBL ID','Smiles']]
    
    smi = df['Smiles']
    sd = [Chem.MolFromSmiles(m) for m in smi]
    y = df['pIC50']
    
        
    desc2d = description2D_calc(sd)
    pd.concat([df_id, y, desc2d],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_2d.csv')),index=False)
    
    estate = estateFP_calc(sd)
    pd.concat([df_id, y, estate],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_estate.csv')),index=False)
    
    autocorr = autocorr2D_calc(sd)
    pd.concat([df_id, y, autocorr],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_autocorr2d.csv')),index=False)
    
    maccs = maccsFP_calc(sd)
    pd.concat([df_id, y, maccs],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_maccs.csv')),index=False)
    
    avalon = avalonFP_calc(sd)
    pd.concat([df_id, y, avalon],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_avalon.csv')),index=False)
    
    avalon_count = avalonCountFP_calc(sd)
    pd.concat([df_id, y, avalon_count],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_avalon_count.csv')),index=False)
    
    layer = layerFP_calc(sd)
    pd.concat([y, layer],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_layer.csv')),index=False)
    
    morgan2 = morganFP_calc(sd, 2)
    pd.concat([y, morgan2],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_morgan2.csv')),index=False)
    
    morgan3 = morganFP_calc(sd, 3)
    pd.concat([y, morgan3],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_morgan3.csv')),index=False)
    
    morgan4 = morganFP_calc(sd, 4)
    pd.concat([y, morgan4],axis=1).to_csv(os.path.join(wdir,t.replace('.csv','_morgan4.csv')),index=False)
    