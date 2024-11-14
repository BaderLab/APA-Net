import torch
import torch.nn as nn
import torch.nn.functional as F


class ConvBlock(nn.Module):
    """
    Convolutional Block for neural networks.
    Args:
        in_channel (int): Number of input channels.
        out_channel (int): Number of output channels.
        ...
    """

    def __init__(
        self,
        in_channel,
        out_channel,
        cnvks=1,
        cnvst=1,
        poolks=1,
        poolst=1,
        pdropout=0,
        activation_t="none",
    ):
        super(ConvBlock, self).__init__()
        activations = {
            "ELU": nn.ELU(),
            "LeakyReLU": nn.LeakyReLU(),
            "none": nn.Identity(),
        }
        self.op = nn.Sequential(
            nn.Conv1d(in_channel, out_channel, cnvks, cnvst, padding=cnvks // 2),
            nn.BatchNorm1d(out_channel),
            activations[activation_t],
            nn.MaxPool1d(kernel_size=poolks, stride=poolst),
            nn.Dropout(p=pdropout),
        )

    def forward(self, x):
        return self.op(x)


class FCBlock(nn.Module):
    """
    Fully Connected Block for neural networks.
    Args:
        layer_dims (list): Dimensions of layers.
        dropouts (list): Dropout values for layers.
        dropout (bool): Whether to apply dropout.
    """

    def __init__(self, layer_dims, dropouts, dropout=False):
        super(FCBlock, self).__init__()
        layers = []
        for i in range(len(layer_dims) - 1):
            layers.append(nn.Linear(layer_dims[i], layer_dims[i + 1]))
            if i < len(layer_dims) - 2:
                layers.append(nn.BatchNorm1d(num_features=layer_dims[i + 1]))
                layers.append(nn.ReLU())
                if dropout:
                    layers.append(nn.Dropout(p=dropouts[i]))
        self.op = nn.Sequential(*layers)

    def forward(self, x):
        return self.op(x)


class ProcessSelfAttn(nn.Module):
    """
    Implements the self-attention mechanism.
    Attributes
    ----------
    nhead : int
        The number of attention heads.
        Each head computes a separate attention score for each token.
    """

    def __init__(
        self,
        embed_dim: int,
        num_layers: int,
        nhead: int,
        dim_feedforward: int = 2048,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.encoder_layer = nn.TransformerEncoderLayer(
            embed_dim,
            nhead,
            dim_feedforward,
            dropout,
            activation="gelu",
            batch_first=True,
        )
        self.transformer = nn.TransformerEncoder(self.encoder_layer, num_layers)

    def forward(self, latent):
        return self.transformer(latent)