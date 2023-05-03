## Train

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import os
import torch
import numpy as np
from torchmetrics.functional import mean_squared_error

## Device Setting
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

## Hyperparameter
SAVE_DIR = "$PATH_TO_SAVE_MODEL"
BATCH_SIZE = 128
LR = 0.001
NUM_EPOCH = 100
N_LAYERS = 5
N_FEATURES = 256


## Train model
if not os.path.exists(SAVE_DIR):
    os.mkdir(SAVE_DIR)

# Set model
from GCN import GCNModel
from Dataset import GCN_data_loaders
model = GCNModel(n_feature=N_FEATURES, num_conv_layers=N_LAYERS)

# Set optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=LR)

import time, copy
train_rmse_history, valid_rmse_history = [], []

best_rmse = 1e6
for i in range(1, NUM_EPOCH + 1):
    ince = time.time()
    # Training phase
    model.train()
    train_batch_rmse = []
    data_loader = GCN_data_loaders['train']

    for batch_idx, batch in enumerate(data_loader):
        x_batch = batch["input"].float().to(device)
        y_batch = batch["output"].float().to(device)
        l_batch = batch["adj"].float().to(device)

        y_pred = model(x_batch, l_batch)
        rmse = torch.sqrt(mean_squared_error(y_pred, y_batch))
        train_batch_rmse.append(rmse.data.cpu().numpy())


        optimizer.zero_grad()
        rmse.backward()
        optimizer.step()

    model.eval()
    with torch.no_grad():
        valid_batch_rmse = []
        data_loader = GCN_data_loaders['val']

        for batch_idx, batch in enumerate(data_loader):
            x_batch = batch["input"].float().to(device)
            y_batch = batch["output"].float().to(device)
            l_batch = batch["adj"].float().to(device)

            y_pred = model(x_batch, l_batch)
            rmse = torch.sqrt(mean_squared_error(y_pred, y_batch))
 
    train_avg_rmse = np.mean(np.array(train_batch_rmse))
    valid_avg_rmse = np.mean(np.array(valid_batch_rmse))

    train_rmse_history.append(train_avg_rmse)
    valid_rmse_history.append(valid_avg_rmse)

    if valid_avg_rmse < best_loss:
            best_epoch = i
            best_loss = valid_avg_rmse

    print(f"\t{i}th EPOCH --- TRAIN RMSE: {train_avg_rmse:.3f} || VALID RMSE: {valid_avg_rmse:.3f} || BEST EPOCH: {best_epoch}", flush=True)

    torch.save(model.state_dict(), os.path.join(SAVE_DIR, f"save_{i}.pt"))
