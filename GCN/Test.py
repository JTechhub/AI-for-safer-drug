## Test

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import os
import torch
from GCN import GCNModel
from Dataset import clean_data
from Dataset import clean_label
from config import df2, BATCH_SIZE, N_FEATURES, N_LAYERS, device
from Dataset import GCNDataset
from torch.utils.data import DataLoader
from Collate_function import sample_collate_fn
from torchmetrics.functional import r2_score
from torchmetrics.functional import mean_absolute_error
from torchmetrics.functional import mean_squared_error
from torchmetrics.functional import spearman_corrcoef
import torchmetrics.functional as Fm

# Set the path
best_model_save_path = "./$PATH_TO_SAVE_MODEL/save_{i}.pt"

if not os.path.exists(best_model_save_path):
    os.mkdir(best_model_save_path)

test_smiles = clean_data(df2)
test_labels = clean_label(df2)

test_loss_history = []

# load dataset
GCN_test_dataset = GCNDataset(test_smiles, test_labels)

# set dataloader
GCN_test_dataloader = DataLoader(GCN_test_dataset, batch_size=BATCH_SIZE, shuffle=False, collate_fn=sample_collate_fn)

save_state_dict = torch.load(best_model_save_path)
model = GCNModel(n_feature=N_FEATURES, num_conv_layers=N_LAYERS).cuda()
model.load_state_dict(save_state_dict)

with torch.no_grad():
    test_batch_losses = []
    for batch_idx, batch in enumerate(GCN_test_dataloader):
        x_batch = batch["input"].float().to(device)
        y_batch = batch["output"].float().to(device)
        l_batch = batch["adj"].float().to(device)

        y_pred = model(x_batch, l_batch)

        loss = torch.sqrt(mean_squared_error(y_pred, y_batch))
        r2 = r2_score(y_pred, y_batch)
        mae = mean_absolute_error(y_batch, y_pred)
        sp_r = Fm.pearson_corrcoef(y_pred, y_batch) if y_pred.shape[1] > 1 else Fm.pearson_corrcoef(y_pred[:, 0], y_batch[:, 0])


print(f"\ TEST R2: {r2:.3f} || TEST LOSS: {loss:.3f} ||TEST MAE: {mae:.3f} || TEST spr: {sp_r:.3f} ", flush=True)