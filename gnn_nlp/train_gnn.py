import torch_geometric
from torch_geometric import data
import os.path as osp

import pandas as pd
import torch
import torch.nn as nn
import numpy as np
from torch_geometric.loader import DataLoader
from torch_geometric.datasets import ZINC, QM9
import random
from tqdm import tqdm

from torchmetrics.functional import r2_score
from torchmetrics.functional import mean_absolute_error
from torchmetrics.functional import mean_squared_error
from torchmetrics.functional import spearman_corrcoef
from utils import *
from models import gnn_model
import matplotlib.pyplot as plt

set_seed(17)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
train = torch.load('data/train_graph_add_del_h')
test = torch.load('data/test_graph_add_del_h')

random.shuffle(train)
train, valid = train[:6000], train[6000:]

batch_size = 128
train_loader = DataLoader(train, batch_size = batch_size, shuffle = True)
valid_loader = DataLoader(valid, batch_size = batch_size, shuffle = False)
test_loader = DataLoader(test, batch_size = batch_size, shuffle = False)

# dataset = QM9('data/')

# pretrain_loader = DataLoader(dataset[:120000], batch_size = 256, shuffle = True)
# prevalid_loader = DataLoader(dataset[120000:], batch_size = 256)

# pretrain_model = gnn_model().to(device)

# criterion = nn.MSELoss()
# optimizer = torch.optim.Adam(pretrain_model.parameters(), lr=1e-3)
# scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10, eta_min=1e-9)

# best_loss = 10e10000000
# train_history = []
# val_history = []
# for epoch in tqdm(range(300)):
    
#     train_loss = 0
#     val_loss = 0

#     pretrain_model.train()
#     for data in pretrain_loader:
#         data = data.to(device)
#         pred = pretrain_model(data, is_pretrain = True)
#         loss = torch.sqrt(criterion(pred, data.y[:,:1]))
#         train_loss+=loss.item()

#         optimizer.zero_grad()
#         loss.backward()
#         optimizer.step()
    
#     # scheduler.step()

#     pretrain_model.eval()
#     for data in prevalid_loader:
#         data = data.to(device)
#         pred = pretrain_model(data, is_pretrain = True)
#         loss = torch.sqrt(criterion(pred, data.y[:,:1]))

#         val_loss+=loss.item()

    

#     train_loss = train_loss/len(train_loader)
#     val_loss = val_loss/len(valid_loader)

#     if val_loss < best_loss:
#         torch.save(pretrain_model.state_dict(), 'model/best_pretrain.pt')
#         best_loss = val_loss

#     train_history.append(train_loss)
#     val_history.append(val_loss)

#     print(f"Epoch {epoch+1} / train loss : {train_loss:.4f}, valid loss : {val_loss:.4f}")

# del pretrain_model, loss

main_model = gnn_model().to(device)

# main_model.load_state_dict(torch.load('model/best_pretrain.pt'))

criterion = nn.MSELoss()
optimizer = torch.optim.Adam(main_model.parameters(), lr=1e-3)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10, eta_min=1e-9)

best_loss = 10e10000000
train_history = []
val_history = []
for epoch in tqdm(range(300)):

    train_loss = 0
    val_loss = 0

    main_model.train()
    for data in train_loader:
        data = data.to(device)
        pred = main_model(data)
        loss = torch.sqrt(criterion(pred, data.y.unsqueeze(1)))
        train_loss+=loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    
    # scheduler.step()

    main_model.eval()
    for data in valid_loader:
        data = data.to(device)
        pred = main_model(data)
        loss = torch.sqrt(criterion(pred, data.y.unsqueeze(1)))

        val_loss+=loss.item()

    

    train_loss = train_loss/len(train_loader)
    val_loss = val_loss/len(valid_loader)

    if val_loss < best_loss:
        torch.save(main_model.state_dict(), 'model/best_gnn_model.pt')
        best_loss = val_loss

    train_history.append(train_loss)
    val_history.append(val_loss)

    print(f"Epoch {epoch+1} / train loss : {train_loss:.4f}, valid loss : {val_loss:.4f}")

del main_model, loss

rmse = 0
r2 = 0
mae = 0
sp_r = 0
best_model = gnn_model().to(device)
best_model.load_state_dict(torch.load('model/best_gnn_model.pt'))
best_model.eval()
for data in test_loader:
    data = data.to(device)
    pred = best_model(data)

    rmse += torch.sqrt(mean_squared_error(pred, data.y.unsqueeze(1))).item()
    r2 += r2_score(pred, data.y.unsqueeze(1)).item()
    mae += mean_absolute_error(pred, data.y.unsqueeze(1)).item()
    sp_r += spearman_corrcoef(pred.squeeze(1), data.y).item()
    
rmse = rmse/len(test_loader)
r2 = r2/len(test_loader)
mae = mae/len(test_loader)
sp_r = sp_r/len(test_loader)

data = {'rmse' : [rmse],
        'r2' : [r2],
        'mae' : [mae],
        'spear' : [sp_r],
        }

df = pd.DataFrame(data)

df.to_csv("gnn_result.csv")