## GCN Code ##

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import os
import torch
import torch.nn as nn
import torch.nn.functional as F

class GCNModel(nn.Module):

    def __init__(self, n_feature, num_conv_layers):
        super(GCNModel, self).__init__()
        self.n_feature = n_feature
        self.num_conv_layers = num_conv_layers
        self.activation = nn.ReLU()

        layer_list = []
        for i in range(num_conv_layers):
          layer_list.append(nn.Linear(n_feature, n_feature))
        self.embedding = nn.Linear(27, n_feature)
        self.gnn = nn.ModuleList([nn.Linear(n_feature, n_feature) for _ in \
                range(self.num_conv_layers)])
        self.readout = nn.ModuleList(layer_list)
        self.fc_layer = nn.Linear(n_feature, 1)


    def forward(self, x, A):
        x = self.embedding(x)

        for l in self.readout:
          x = l(x)     
          x = torch.matmul(A,x)
          x = F.relu(x)
        x = x.mean(1)
        retval = self.fc_layer(x)
        return retval  
    
