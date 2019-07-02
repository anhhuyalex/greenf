import torch
import torch.nn as nn
import torch.autograd as autograd
from torch.autograd import Variable
import torch.nn.functional as F
import math
import numpy as np
import sys; 
from tensorboardX import SummaryWriter
import torchsummary
sys.path


# #parameters
num_res_layers = 7
cnn_num_filters = 16
cnn_kernel_size = 5
cnn_stride_size = 1
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')



def form_hamiltonian(sites=5, e_n = 1, t_n = 1):
    
    def gen_triples(N):
        i, j, k = (0, 0, 0)
        coords = []
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    coords.append((i,j,k))

        return coords
    
    coords = gen_triples(sites)
    mapping = {i: coords[i] for i in range(len(coords))}
    inv_mapping = {coords[i]: i for i in range(len(coords))}
    H = np.zeros((sites ** 3, sites ** 3), dtype = complex)
    print(len(coords))
    for i in range(sites ** 3):
        x, y, z = mapping[i]
        H[i, i] = e_n
        if (z + 1 < sites):
            H[i, i+1] = t_n
        if y + 1 < sites:
            H[i, i+sites] = t_n
        if x + 1 < sites:
            H[i, i+sites**2] = t_n
        if z > 0:
            H[i, i-1] = t_n
        if y > 0:
            H[i, i-sites] = t_n
        if x > 0:
            H[i, i-sites**2] = t_n
            
    return gen_triples, H

gen_triples, H = form_hamiltonian()
H

H = torch.Tensor(np.real(H))
z = torch.zeros_like(H).fill_(3)
inp = torch.cat((H.unsqueeze(0), z.unsqueeze(0)), dim = 0)
inp = inp.unsqueeze(0)
inp.shape

from collections import OrderedDict
class CNNEncoder(nn.Module):
    """Maps a Hamiltonian and an energy input a hidden vector"""
    def __init__(self):
        super(CNNEncoder, self).__init__()
        self.layer1 = nn.Sequential(
            nn.Conv2d(2, cnn_num_filters, kernel_size=cnn_kernel_size, stride=cnn_stride_size, padding=2),
            nn.BatchNorm2d(cnn_num_filters),
            nn.ReLU())
        self.test = nn.Sequential(OrderedDict([
            ("first_conv", nn.Conv2d(2, cnn_num_filters, kernel_size=cnn_kernel_size, stride=cnn_stride_size, padding=2)),
            ("first_bn", nn.BatchNorm2d(cnn_num_filters))
        ]))
       
        self.resblock = [self._build_res_block(index) for index in range(cnn_num_filters)]
        self.final_encoder = nn.Sequential(OrderedDict([
            ("last_conv", nn.Conv2d(cnn_num_filters, 2, kernel_size=1, stride=cnn_stride_size, padding=2)),
            ("last_bn", nn.BatchNorm2d(2)),
            ("last_relu", nn.ReLU())
        ]))

    def forward(self, inputs):
        out = self.layer1(inputs)
        for _ in range(num_res_layers):
          res = out.clone()
          resblock = self.resblock[_]
          out = resblock(out) + res
        out = self.final_encoder(out)
        return out
      
    def _build_res_block(self, index):
        resblock = nn.Sequential(OrderedDict([
            ("conv1_resblock"+str(index), nn.Conv2d(cnn_num_filters, cnn_num_filters, kernel_size=cnn_kernel_size, stride=cnn_stride_size, padding=2)),
            ("bn1_resblock"+str(index), nn.BatchNorm2d(cnn_num_filters)),
            ("relu1_resblock"+str(index), nn.ReLU()),
            ("conv2_resblock"+str(index), nn.Conv2d(cnn_num_filters, cnn_num_filters, kernel_size=cnn_kernel_size, stride=cnn_stride_size, padding=2)),
            ("bn2_resblock"+str(index), nn.BatchNorm2d(cnn_num_filters)),
            ("relu1_resblock"+str(index), nn.ReLU())]))
        return resblock
    
v = CNNEncoder().to(device)
torchsummary.summary(v, (2, 125, 125))
        
        
        
    
    
w = SummaryWriter(log_dir="./logs", comment='CNNEncoder') 
w.add_graph(v, inp)
w.close()
