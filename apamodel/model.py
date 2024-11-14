import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from blocks import ConvBlock, FCBlock, ProcessSelfAttn
import numpy as np

RBP_COUNT = 327
FIX_SEQ_LEN = 4000



class APAData(Dataset):
    """
    APAData is a dataset class for APA-Net model.
    Args:
        seqs (Tensor): Sequences tensor.
        df (DataFrame): Dataframe containing sample information.
        ct (DataFrame): Cell type profiles.
        device (str): Device to use (e.g., 'cuda' or 'cpu').
    """

    def __init__(self, data, device):
        self.device = device
        self.reg_label = torch.from_numpy(
            np.array(data[:, 7].tolist(), dtype=np.float32)
        ).to(device)
        self.oneH_seqs = torch.from_numpy(np.array(data[:, 6].tolist())).to(
            device
        )
        self.ct_profiles = torch.from_numpy(np.array(data[:, 8].tolist())).to(
            device
        )
        self.celltype_name = data[:, 1].tolist()
        self.switch_name = data[:, 5].tolist()

    def __len__(self):
        return self.reg_label.shape[0]

    def __getitem__(self, idx):
        seq = self.oneH_seqs[idx].type(torch.cuda.FloatTensor)
        reg_label = self.reg_label[idx]
        celltype_profile = self.ct_profiles[idx].type(torch.cuda.FloatTensor)
        celltype_name = self.celltype_name[idx]
        switch_name = self.switch_name[idx]
        return (seq, reg_label, celltype_profile, celltype_name, switch_name)

class APANET(nn.Module):
    """
    APANET is a deep neural network for APA-Net.
    Includes Convolutional, Attention, and Fully Connected blocks.
    """

    def __init__(self, config):
        super(APANET, self).__init__()
        self.config = config
        self.device = config["device"]
        self._build_model()

    def _build_model(self):
        # Convolutional Block
        self.conv_block_1 = ConvBlock(
            in_channel=4,
            out_channel=self.config["conv1kc"],
            cnvks=self.config["conv1ks"],
            cnvst=self.config["conv1st"],
            poolks=self.config["pool1ks"],
            poolst=self.config["pool1st"],
            pdropout=self.config["cnvpdrop1"],
            activation_t="ELU",
        )
        # Calculate output length after Convolution
        cnv1_len = self._get_conv1d_out_length(
            FIX_SEQ_LEN,
            self.config["conv1ks"],
            self.config["conv1st"],
            self.config["pool1ks"],
            self.config["pool1st"],
        )

        # Attention Block
        self.attention = nn.MultiheadAttention(
            embed_dim=self.config["conv1kc"],
            num_heads=self.config["Matt_heads"],
            dropout=self.config["Matt_drop"],
        )

        # Fully Connected Blocks
        fc1_L1 = cnv1_len * self.config["conv1kc"]
        self.fc1 = FCBlock(
            layer_dims=[fc1_L1, *self.config["fc1_dims"]],
            dropouts=self.config["fc1_dropouts"],
            dropout=True,
        )

        fc2_L1 = self.config["fc1_dims"][-1] + RBP_COUNT
        self.fc2 = FCBlock(
            layer_dims=[fc2_L1, *self.config["fc2_dims"]],
            dropouts=self.config["fc2_dropouts"],
            dropout=True,
        )
        self.process_self_attn = ProcessSelfAttn(
            self.config["psa_query_dim"],
            self.config["psa_num_layers"],
            self.config["psa_nhead"],
            self.config["psa_dim_feedforward"],
            self.config["psa_dropout"]
        )

    def _get_conv1d_out_length(self, l_in, kernel, stride, pool_kernel, pool_stride):
        """Utility method to calculate output length of Conv1D layer."""
        length_after_conv = (
            l_in + 2 * (kernel // 2) - 1 * (kernel - 1) - 1
        ) // stride + 1
        return (length_after_conv - pool_kernel) // pool_stride + 1

    def forward(self, seq, celltype):
        # Convolutional forward
        x_conv = self.conv_block_1(seq) # batch, 64/128(dim), 80(len)
        x = x_conv.permute(0, 2, 1)  # reshape for attention block so dim is first
        # x, _ = self.attention(x, x, x)
        x = self.process_self_attn(x)
        x = x.permute(0, 2, 1)  # reshape back
        x = x + x_conv  # add residual connection
        x = torch.flatten(x, 1)  # flatten for FC layers
        x = self.fc1(x)  # FC block 1
        x = torch.cat((x, celltype), 1)  # concat with celltype profile
        x = self.fc2(x)  # FC block 2
        return x

    def compile(self):
        """Compile the model with optimizer and loss function."""
        self.to(self.device)
        if self.config["opt"] == "Adam":
            self.optimizer = optim.AdamW(
                self.parameters(),
                weight_decay=self.config["adam_weight_decay"],
                lr=self.config["lr"],
            )
        if self.config["loss"] == "mse":
            self.loss_fn = nn.MSELoss()

    def save_model(self, filename):
        torch.save(self.state_dict(), filename)

    def load_model(self, filename):
        self.load_state_dict(torch.load(filename))
