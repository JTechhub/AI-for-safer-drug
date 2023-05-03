## Preprocessing

"""
Created on Wed May 3 2023

@author: Minju Na
"""

from config import df1

# Preprocessing
def clean_id(data):
    id_list = []
    for i in range(len(data)):
        id_list.append(data.iloc[i,0])
    return id_list


def clean_data(data):
    new_smiles_list = []
    for i in range(len(data)):
      if len(data.iloc[i,1]) > 0:
          new_smiles_list.append(data.iloc[i,1])
    return new_smiles_list


def clean_label(data):
    label_list = []
    for i in range(len(data)):
        label_list.append(data.iloc[i,2])
    return label_list


id = clean_id(df1)
smiles_list = clean_data(df1)
label_list = clean_label(df1)