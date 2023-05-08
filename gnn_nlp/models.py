
from torch_geometric.nn import DimeNet, DimeNetPlusPlus
from torch_geometric.nn.norm import BatchNorm
import torch.nn as nn

class gnn_model(nn.Module):
    def __init__(self):
        super().__init__()
        self.dime = DimeNetPlusPlus(hidden_channels=128, out_channels=16, num_blocks=4, int_emb_size=64,
                      basis_emb_size=8, out_emb_channels=256, num_spherical=7, num_radial=6)
        self.lr = nn.LeakyReLU()
        self.drop = nn.Dropout(p=0.2)
        
        self.bn = BatchNorm(16)
        self.linear = nn.Linear(16,1)

        self.bn_pretrain = BatchNorm(16)
        self.pretrain = nn.Linear(16,1)

    def forward(self, data, is_pretrain = False):
        if is_pretrain:
            x = self.dime(data.z,data.pos,data.batch)
            x = self.bn_pretrain(x)
            x = self.lr(x)
            x = self.drop(x)
            x = self.pretrain(x)
        else:
            x = self.dime(data.x[:,0],data.pos,data.batch)
            x = self.bn(x)
            x = self.lr(x)
            x = self.drop(x)
            x = self.linear(x)
        return x
    
