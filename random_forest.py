import sys
 
import os
import rdkit

import pandas as pd
import numpy as np

from rdkit_descriptor_calculation import *
from scipy import stats

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split

from scipy.stats import spearmanr
from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import r2_score
import math

train_df = pd.read_csv("../data/hERG_train.csv")
test_df = pd.read_csv("../data/hERG_test.csv")

target = [train_df, test_df]

train_set = []
test_set = []

sets = [train_set, test_set]

targetStr = ["train","test"]

for i, df in enumerate(target):
    df_id = df[['Molecule ChEMBL ID','Smiles']]
    
    smi = df['Smiles']
    sd = [Chem.MolFromSmiles(m) for m in smi]
    y = df['pIC50']
    
    print(1)
    desc2d = pd.read_csv("../data/hERG_"+ targetStr[i] + "_2d.csv")            
    sets[i].append(desc2d)
    print(2)
    
    estate = pd.read_csv("../data/hERG_"+ targetStr[i] + "_estate.csv")    
    sets[i].append(estate)
    print(3)        
    
    autocorr = pd.read_csv("../data/hERG_"+ targetStr[i] + "_autocorr2d.csv")    
    sets[i].append(autocorr)
    print(4)        
    
    maccs = pd.read_csv("../data/hERG_"+ targetStr[i] + "_maccs.csv")    
    sets[i].append(maccs)

    print(5)        
    avalon = pd.read_csv("../data/hERG_"+ targetStr[i] + "_avalon.csv")    
    sets[i].append(avalon)
    
    print(6)        
    avalon_count = pd.read_csv("../data/hERG_"+ targetStr[i] + "_avalon_count.csv")    
    sets[i].append(avalon_count)
    
    print(7)        
    layer = pd.read_csv("../data/hERG_"+ targetStr[i] + "_layer.csv")    
    sets[i].append(pd.concat([df_id, layer],axis=1))    
    
    print(8)        
    morgan2 = pd.read_csv("../data/hERG_"+ targetStr[i] + "_morgan2.csv")    
    sets[i].append(pd.concat([df_id, morgan2],axis=1))    
    
    print(9)       
    morgan3 = pd.read_csv("../data/hERG_"+ targetStr[i] + "_morgan3.csv")    
    sets[i].append(pd.concat([df_id, morgan3],axis=1))
    
    print(10)        
    morgan4 = pd.read_csv("../data/hERG_"+ targetStr[i] + "_morgan4.csv")        
    sets[i].append(pd.concat([df_id, morgan4],axis=1))    
    
    print(11)        

desc = ["desc2d","estateFP","autocorr2D","maccsFP","avalonFP","avalonCountFP","layerFP","morgan2","morgan3","morgan4"]
     
len(desc)    
    
# Fitting the model
for i in range(len(desc)):
    X_train = train_set[i].drop(['Molecule ChEMBL ID','Smiles', 'pIC50'], axis=1)
    y_train = train_set[i]['pIC50']

    X_test = test_set[i].drop(['Molecule ChEMBL ID','Smiles', 'pIC50'], axis=1)
    y_test = test_set[i]['pIC50']

    print("* Descriptor Name : %s" %(desc[i]))
    print("  Training...")

    # find best hyper params
    
    params = {    
      'n_estimators':(100, 200),
      'max_depth' : (5, 8),
      'min_samples_leaf' : (8, 18),
      'min_samples_split' : (8, 16)            
      }
    
    rf_run = RandomForestRegressor(random_state=0, n_jobs=-1)
    grid_cv = GridSearchCV(rf_run, param_grid=params, cv=2, n_jobs=-1)
    grid_cv.fit(X_train, y_train)
 
    print('best hyper params:', grid_cv.best_params_)

    # fit
    rf_run = RandomForestRegressor(random_state=0, max_depth=grid_cv.best_params_['max_depth'], min_samples_leaf=grid_cv.best_params_['min_samples_leaf'], min_samples_split=grid_cv.best_params_['min_samples_split'],n_estimators=grid_cv.best_params_['n_estimators'])
    rf_run.fit(X_train, y_train)
 
    # Predict the model
    pred = rf_run.predict(X_test)
  
    r_sq = r2_score(y_test, pred)
    rmse = np.sqrt(MSE(y_test, pred))
    sp_r = stats.spearmanr(y_test, pred).statistic
  
    print(" > Evaluation Metrics") 
    print("   R square : %f" %(r_sq))
    print("   RMSE : % f" %(rmse))
    print("   Spearman R : %f" %(sp_r))
    print("\n")
    