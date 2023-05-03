## config

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import os
import torch
import pandas as pd

train_csv = "train.csv"
test_csv = "test.csv"

# Set smiles & label
df1 = pd.read_csv(train_csv)
index = df1["index"].to_numpy()
smiles = df1["Smiles"].to_numpy()
labels = df1["pIC50"].to_numpy()

df2 = pd.read_csv(test_csv)
test_index = df2["index"].to_numpy()
test_smiles = df2["Smiles"].to_numpy()
test_labels = df2["pIC50"].to_numpy()

## Hyperparameter
SAVE_DIR = "$PATH_TO_SAVE_MODEL"
BATCH_SIZE = 128
LR = 0.001
NUM_EPOCH = 100
N_LAYERS = 5
N_FEATURES = 256

# Settings
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

os.chdir(f'/Set you dir/')
