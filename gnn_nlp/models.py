
from torch_geometric.nn import DimeNet, DimeNetPlusPlus
import torch.nn as nn
from transformers import  AutoModelForMaskedLM

class gnn_model(nn.Module):
    def __init__(self):
        super().__init__()
        self.dime = DimeNetPlusPlus(hidden_channels=128, out_channels=16, num_blocks=4, int_emb_size=64,
                      basis_emb_size=8, out_emb_channels=256, num_spherical=7, num_radial=6)
        self.lr = nn.LeakyReLU()
        self.drop = nn.Dropout(p=0.5)
        
        self.linear = nn.Linear(16,1)

    def forward(self, data):
            x = self.dime(data.x[:,0],data.pos,data.batch)
            x = self.lr(x)
            x = self.drop(x)
            x = self.linear(x)
            return x
    


class nlp_model(nn.Module):
    def __init__(self):
        super().__init__()
        self.roberta = AutoModelForMaskedLM.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k")

        self.lr = nn.LeakyReLU()
        self.drop = nn.Dropout(p=0.5)

        self.linear1 = nn.Linear(52000,1024)
        self.linear2 = nn.Linear(1024,32)
        self.linear3 = nn.Linear(32,1)

    def forward(self, d0, d1):
        x = self.roberta(d0,d1)
        x = x.logits[:,0,:]

        x = self.linear1(x)
        x = self.lr(x)
        x = self.drop(x)
        x = self.linear2(x)
        x = self.lr(x)
        x = self.drop(x)
        x = self.linear3(x)
        return x