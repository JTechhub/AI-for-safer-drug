## Define Dataset

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import torch
import config
import numpy as np
from functools import reduce
from rdkit import Chem
from rdkit.Chem import AllChem

class GCNDataset(Dataset):
  
  def __init__(self, smiles_list, label_list):
    super().__init__()
    self.smiles_list = smiles_list
    self.label_list = label_list
    self.id_list = id
    
    self._set_mol_list()
    self._set_atom_list()

    self.atom_feature = [] 
    self.adj_list = []

  def __len__(self):
    return len(self.smiles_list)

  def __getitem__(self,idx):
    sample = dict()
    mol = self.mol_list[idx]
    l = self.label_list[idx]
    if self.label_list is not None:
        label = self.label_list[idx]
    else:
        label = 0.

    input = self._get_node_feature_matrix(mol)
    output = torch.tensor([l])
    adj = self._get_adjacency_matrix(mol)

    sample = {
                "input": torch.Tensor(input), 
                "output": torch.Tensor(output),
                "adj": torch.Tensor(adj), 
                "num_atom": torch.LongTensor([mol.GetNumAtoms()])
    }
    return sample

  def _set_mol_list(self):
        assert len(self.smiles_list) > 0
        self.mol_list = [Chem.MolFromSmiles(smi) for smi in self.smiles_list]

  def _set_atom_list(self, atom_list=None):
    assert len(self.mol_list) > 0
    if atom_list is not None:
        self.atom_list = atom_list
    else:
        whole_atom = list(reduce(lambda x, y: x | y, \
                                 [set([a.GetSymbol() for a in mol.GetAtoms()]) \
                                 for mol in self.mol_list]))
        self.atom_list = sorted(whole_atom)

  def _get_num_atom_feature(self):
      assert len(self.atom_list) > 0
      return len(self.atom_list)

  def _get_node_feature_matrix(self, mol):
      return np.array([self._get_atom_feature_vector(a) for a in mol.GetAtoms()])

  def _get_adjacency_matrix(self, mol):
        """
        return normalized adjacency matrix of a given molecule
        use GetAdjacencyMatrix()
        """
        def normalizeAdjacency(adj):
            """
            return a normalized adjacency matrix: D^-1/2 @ (A + I) @ D^-1/2
            """
            assert adj.shape[0] == adj.shape[1]
            A = adj + np.eye(adj.shape[0]) 
            d = np.sum(A, axis=1) 
            d = 1/np.sqrt(d)
            D = np.diag(d)
            return D @ A @ D

        adj = AllChem.GetAdjacencyMatrix(mol)
        return normalizeAdjacency(adj)

  def one_of_k_encoding(self, x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

  def _get_atom_feature_vector(self, atom): 
      """
      return one-hot encoded vector of a given atom
      use atom.GetSymbol() and self._get_one_hot_vector()
      """
      return np.array(self.one_of_k_encoding(atom.GetSymbol(),['C','N','O','F','S','P','Cl','Br','ELSE']) +
                      self.one_of_k_encoding(atom.GetTotalNumHs(), [0, 1, 2, 3, 4,'ELSE']) +
                      self.one_of_k_encoding(atom.GetExplicitValence(), [0, 1, 2, 3, 4,'ELSE']) +
                      self.one_of_k_encoding(atom.GetFormalCharge(),[-2, -1,0,1,2,'ELSE']))


  def _get_one_hot_vector(self, idx, vec_dim):
      return np.eye(vec_dim)[idx]
  

# Split train / valid
from sklearn.model_selection import train_test_split
import numpy as np
random_seed = 5
np.random.seed(random_seed)
train_smiles, valid_smiles, train_labels, valid_labels = train_test_split(config.smiles, config.labels, test_size=0.2, shuffle=True)
print(f"train_smiles_list: {len(train_smiles)}, valid_smiles_list: {len(valid_smiles)}")


# Dataloader
# GCN dataset
GCN_train_data = GCNDataset(train_smiles, train_labels)
GCN_val_data = GCNDataset(valid_smiles, valid_labels)

from torch.utils.data import DataLoader
from Collate_function import sample_collate_fn
GCN_data_loaders = {}
GCN_data_loaders['train'] = DataLoader(GCN_train_data, batch_size = config.BATCH_SIZE, shuffle = True, collate_fn=sample_collate_fn)
GCN_data_loaders['val'] = DataLoader(GCN_val_data, batch_size = config.BATCH_SIZE,shuffle = False, collate_fn=sample_collate_fn)
