"""
This module contains implementations of the ML models that were presented in 
the study of disorder in computational (materials) databases.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class LinDisorderClassifier(nn.Module):
    def __init__(self, nin, nout, nh) -> None:
        super().__init__()
        self.fc1 = nn.Linear(nin, nh)
        self.fc2 = nn.Linear(nh, nout)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.sigmoid(self.fc2(x))
        return x


class RNNDisorderClassifier(nn.Module):
    def __init__(self, nin, nh, nl, nout, batched=False) -> None:
        super().__init__()
        # Dimensions
        self.input_dim = nin
        self.hidden_dim = nh
        self.layer_dim = nl
        self.output_dim = nout
        # Batch size 1?
        self.batched = batched
        # Central LSTM
        if self.batched:
            self.lstm = torch.nn.LSTM(nin, nh, nl, batch_first=True)
        else:
            self.lstm = torch.nn.LSTM(nin, nh, nl)
        self.lstm.double()
        # Readout layer
        self.fc1 = torch.nn.Linear(nh, nh)
        self.fc2 = torch.nn.Linear(nh, nout)
        self.fc1.double()
        self.fc2.double()
    
    def forward(self, x) -> torch.Tensor:
        if self.batched:
            h0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim).requires_grad_()
            c0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim).requires_grad_()
        else:
            h0 = torch.zeros(self.layer_dim, self.hidden_dim).requires_grad_()
            c0 = torch.zeros(self.layer_dim, self.hidden_dim).requires_grad_()
        out, (hn, cn) = self.lstm(x.double(), (h0.double(), c0.double()))
        out = F.relu(self.fc1(out.double()))
        out = F.sigmoid(self.fc2(out.double()))[:,-1]
        return out


class RNNDisorderClassifier_general(nn.Module):
    def __init__(self, nin, nh, nl, nout, nread=1, batched=False) -> None:
        super().__init__()
        # Dimensions
        self.input_dim = nin
        self.hidden_dim = nh
        self.layer_dim = nl
        self.output_dim = nout
        self.num_readout = nread
        # Batch size 1?
        self.batched = batched
        # Central LSTM
        if self.batched:
            self.lstm = torch.nn.LSTM(nin, nh, nl, batch_first=True)
        else:
            self.lstm = torch.nn.LSTM(nin, nh, nl)
        self.lstm.double()
        # Readout layer
        rls = [
            torch.nn.Linear(nh,nh) for _ in range(self.num_readout - 1)
        ]
        layers = []
        for l in rls:
            layers.append(l)
            layers.append(torch.nn.ReLU())
        layers.append(torch.nn.Linear(nh,nout))
        self.fc = torch.nn.Sequential(*layers)
        self.fc.double()
    
    def forward(self, x) -> torch.Tensor:
        if self.batched:
            h0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim).requires_grad_()
            c0 = torch.zeros(self.layer_dim, x.size(0), self.hidden_dim).requires_grad_()
        else:
            h0 = torch.zeros(self.layer_dim, self.hidden_dim).requires_grad_()
            c0 = torch.zeros(self.layer_dim, self.hidden_dim).requires_grad_()
        out, (hn, cn) = self.lstm(x.double(), (h0.double(), c0.double()))
        out = F.sigmoid(self.fc(out.double()))[:,-1]
        return out
    

### NOTE: this is a deterministic model that has NO learnable parameters!
class BaselineClassifier(nn.Module):
    POOLINGS = {
        'max': torch.amax,
        'softmax': torch.softmax
    }
    def __init__(self, pairprobs, pool='max'):
        super().__init__()
        # Store mapping elpair -> probability
        if isinstance(pairprobs, torch.Tensor):
            self.pairprobs = pairprobs          # Either as tensor input
        else:
            torch.load(pairprobs)               # Or read from file
        self.pool = self.POOLINGS[pool]
    
    def forward(self, x) -> torch.Tensor:
        # 1. Assign probability to each element pair
        # -> take 2D matrix rep & multiply it with our elpair matrix
        x = torch.nan_to_num(
            self.pairprobs[torch.newaxis,:,:] * x,
            nan=0.
        )
        # 2. Apply pooling
        # -> max pooling
        return self.pool(x, dim=(1, 2))