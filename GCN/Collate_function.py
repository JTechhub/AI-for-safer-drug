## Collate function

"""
Created on Wed May 3 2023

@author: Minju Na
"""

import torch

def sample_collate_fn(samples):
    '''
    Dataloader will make a list of samples with a len(samples) = batch_size.
    Collate function should pad all the tensors in every sample at maximum size,
    and stack them on a batch dimension.
    '''
    def get_max_size(tensors):
        size_list = [tensor.shape for tensor in tensors]
        max_size = sorted(size_list, \
                key=lambda x: int(torch.prod(torch.tensor(x))))[-1]
        return max_size

    def max_padding(tensor, size):
        assert tensor.shape[0] <= size[0] and tensor.shape[1] <= size[1]
        max_tensor = torch.zeros(*size)
        max_tensor[:tensor.shape[0], :tensor.shape[1]] = tensor
        return max_tensor

    input_list = [sample["input"] for sample in samples]
    adj_list = [sample["adj"] for sample in samples]
    input_max_size = get_max_size(input_list)
    adj_max_size = get_max_size(adj_list)

    inputs = torch.stack([max_padding(input, input_max_size) \
            for input in input_list], dim=0)
    adjs = torch.stack([max_padding(adj, adj_max_size) \
            for adj in adj_list], dim=0)
    outputs = torch.stack([sample["output"] for sample in samples], dim=0)
    num_atoms = torch.stack([sample["num_atom"] for sample in samples], dim=0)

    sample_batch = {
            "input": inputs,
            "output": outputs,
            "adj": adjs,
            "num_atom": num_atoms
    }
    return sample_batch