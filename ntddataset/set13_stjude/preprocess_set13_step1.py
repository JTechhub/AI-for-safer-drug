# -*- coding: utf-8 -*-
import os
import sys
import string
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
import numpy as np


originalfile = './Identification_of_Selective_Inhibitors_of_PfHT.sdf'
newfile = 'dataset13.sdf'
datafile = open(newfile,'w')

with open(originalfile,'r',errors='replace') as f:
  text = f.read()  
  datafile.write(text)
  f.close()
  datafile.close()

sdf = PandasTools.LoadSDF(newfile)

sdf.to_excel('dataset13.xlsx', # directory and file name to write
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


