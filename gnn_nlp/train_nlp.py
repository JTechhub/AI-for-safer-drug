from transformers import AutoTokenizer
from torch_geometric import data
import os.path as osp

import pandas as pd
import torch
import torch.nn as nn
import numpy as np
import random
from tqdm import tqdm
import copy
from models import nlp_model

from torchmetrics.functional import r2_score
from torchmetrics.functional import mean_absolute_error
from torchmetrics.functional import mean_squared_error
from torchmetrics.functional import spearman_corrcoef

# load pretrained model
tokenizer = AutoTokenizer.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")

train = torch.load('data/train_graph_add_del_h')
test = torch.load('data/test_graph_add_del_h')

random.shuffle(train)
train, valid = train[:6000], train[6000:]

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

nlp_train_x = []
nlp_train_y = []

nlp_valid_x = []
nlp_valid_y = []

nlp_test_x = []
nlp_test_y = []

for d in train:
    nlp_train_x.append(d.smiles)
    nlp_train_y.append(d.y)

for d in valid:
    nlp_valid_x.append(d.smiles)
    nlp_valid_y.append(d.y)

for d in test:
    nlp_test_x.append(d.smiles)
    nlp_test_y.append(d.y)

nlp_train_x = tokenizer(nlp_train_x, padding='max_length', max_length=150, return_tensors="pt")
nlp_train_y = torch.tensor(nlp_train_y).unsqueeze(-1)

nlp_valid_x = tokenizer(nlp_valid_x, padding='max_length', max_length=150, return_tensors="pt")
nlp_valid_y = torch.tensor(nlp_valid_y).unsqueeze(-1)

nlp_test_x = tokenizer(nlp_test_x, padding='max_length', max_length=150, return_tensors="pt")
nlp_test_y = torch.tensor(nlp_test_y).unsqueeze(-1)

train_dataset = torch.utils.data.TensorDataset(nlp_train_x.input_ids, nlp_train_x.attention_mask, nlp_train_y)
valid_dataset = torch.utils.data.TensorDataset(nlp_valid_x.input_ids, nlp_valid_x.attention_mask, nlp_valid_y)
test_dataset = torch.utils.data.TensorDataset(nlp_test_x.input_ids, nlp_test_x.attention_mask, nlp_test_y)

batch_size= 128

train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
valid_loader = torch.utils.data.DataLoader(valid_dataset, batch_size=batch_size)
test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size)

model = nlp_model().to(device)

criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10, eta_min=0)

best_loss = 10e100000
for epoch in tqdm(range(200)):
    train_history = []
    val_history = []
    train_loss = 0
    val_loss = 0

    model.train()
    for data in train_loader:
        data[0], data[1], data[2] = data[0].to(device), data[1].to(device), data[2].to(device)
        pred = model(data[0], data[1])
        
        loss = torch.sqrt(criterion(pred, data[-1]))

        train_loss+=loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    
    scheduler.step()

    model.eval()
    for data in valid_loader:
        data[0], data[1], data[2] = data[0].to(device), data[1].to(device), data[2].to(device)
        pred = model(data[0], data[1])

        loss = torch.sqrt(criterion(pred, data[-1]))

        val_loss+=loss.item()

    train_loss = train_loss/len(train_loader)
    val_loss = val_loss/len(valid_loader)

    if val_loss < best_loss:
        torch.save(model.state_dict(), 'model/best_nlp_model.pt')
        best_loss = val_loss

    train_history.append(train_loss)
    val_history.append(val_loss)

    print(f"Epoch {epoch+1} / train loss : {train_loss}, valid loss : {val_loss}")

del model, loss

rmse = 0
r2 = 0
mae = 0
sp_r = 0
best_model = nlp_model().to(device)
best_model.load_state_dict(torch.load('model/best_nlp_model.pt'))
best_model.eval()
for data in test_loader:
    data[0], data[1], data[2] = data[0].to(device), data[1].to(device), data[2].to(device)
    pred = best_model(data[0], data[1])

    rmse += torch.sqrt(mean_squared_error(pred, data[-1])).item()
    r2 += r2_score(pred, data[-1]).item()
    mae += mean_absolute_error(pred, data[-1]).item()
    sp_r += spearman_corrcoef(pred.squeeze(1), data[-1].squeeze(1)).item()
    
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

df.to_csv("nlp_results/nlp_result.csv")