from tqdm import tqdm
import argparse
import random
import numpy as np
import pandas as pd
import scipy.stats as stats
import torch
from torch.utils.data import DataLoader
from model import APANET, APAData
import wandb


def build_dataloaders(
    device, train_seq, valid_seq, train_data, val_data, batch_size, ct_profiles
):
    """
    Create training and validation data loaders.
    Args:
        device: The device to use for tensors.
        train_seq, valid_seq, train_data, val_data: Data sequences and labels.
        batch_size: Batch size for data loaders.
        ct_profiles: Cell type profiles.
    Returns:
        Tuple of DataLoader for training and validation datasets.
    """
    train_loader = DataLoader(
        APAData(train_seq, train_data, ct_profiles, device),
        batch_size=batch_size,
        shuffle=True,
        drop_last=True,
    )
    valid_loader = DataLoader(
        APAData(valid_seq, val_data, ct_profiles, device),
        batch_size=batch_size,
        shuffle=False,
        drop_last=False,
    )
    return train_loader, valid_loader


def train_one_epoch(model, train_loader):
    """
    Train the model for one epoch.
    Args:
        model: The neural network model.
        train_loader: DataLoader for training data.
    Returns:
        Tuple of average training loss and Pearson correlation coefficient.
    """
    model.train()
    total_loss, predictions, targets = 0.0, [], []
    for seq_X, celltype, _, Y in train_loader:
        model.optimizer.zero_grad()
        outputs = torch.squeeze(model(seq_X, celltype))
        loss = torch.sqrt(model.loss_fn(outputs, Y))
        loss.backward()
        model.optimizer.step()
        total_loss += loss.item() * seq_X.size(0)
        predictions.append(outputs.cpu().detach().numpy())
        targets.append(Y.cpu().detach().numpy())
    avg_loss = total_loss / len(train_loader.dataset)
    correlation = stats.pearsonr(np.concatenate(targets), np.concatenate(predictions))[
        0
    ]
    return avg_loss, correlation


def validate_one_epoch(model, valid_loader):
    """
    Validate the model for one epoch.
    Args:
        model: The neural network model.
        valid_loader: DataLoader for validation data.
    Returns:
        Tuple of average validation loss and Pearson correlation coefficient.
    """
    model.eval()
    total_loss, predictions, targets = 0.0, [], []
    with torch.no_grad():
        for seq_X, celltype, _, Y in valid_loader:
            outputs = torch.squeeze(model(seq_X, celltype))
            loss = torch.sqrt(model.loss_fn(outputs, Y))
            total_loss += loss.item() * seq_X.size(0)
            predictions.append(outputs.cpu().numpy())
            targets.append(Y.cpu().numpy())

    avg_loss = total_loss / len(valid_loader.dataset)
    correlation = stats.pearsonr(np.concatenate(targets), np.concatenate(predictions))[
        0
    ]
    return avg_loss, correlation


def main_train(
    train_seq,
    valid_seq,
    train_data,
    val_data,
    profiles,
    modelfile,
    device,
    config,
    use_wandb,
):
    """
    Main training loop.
    Args:
        train_seq, valid_seq, train_data, val_data: Data sequences and labels.
        profiles: Cell type profiles.
        modelfile: File path to save the best model.
        device: The device to use for tensors.
        config: Configuration parameters for the model.
    """

    use_wandb = args.use_wandb.lower() == "true"
    train_loader, valid_loader = build_dataloaders(
        device,
        train_seq,
        valid_seq,
        train_data,
        val_data,
        config["batch_size"],
        profiles,
    )
    with tqdm(range(config["epochs"]), unit="epoch") as tepochs:
        if use_wandb:
            wandb.login()
            with wandb.init(
                project=config["project_name"],
                settings=wandb.Settings(start_method="thread"),
            ):
                model = APANET(config)
                model.compile()
                best_model_metric = -float("inf")
                for epoch in tepochs:
                    train_loss, train_corr = train_one_epoch(model, train_loader)
                    valid_loss, valid_corr = validate_one_epoch(model, valid_loader)
                    tepochs.set_postfix(train_loss=train_loss, valid_loss=valid_loss)
                    if valid_corr > best_model_metric:
                        best_model_metric = valid_corr
                        model.save_model(modelfile)

                    wandb.log(
                        {
                            "train_loss": train_loss,
                            "train_corr": train_corr,
                            "valid_loss": valid_loss,
                            "valid_corr": valid_corr,
                            "epoch": epoch,
                        }
                    )
        else:
            model = APANET(config)
            model.compile()
            best_model_metric = -float("inf")
            for epoch in tepochs:
                train_loss, train_corr = train_one_epoch(model, train_loader)
                valid_loss, valid_corr = validate_one_epoch(model, valid_loader)
                tepochs.set_postfix(train_loss=train_loss, valid_loss=valid_loss)
                if valid_corr > best_model_metric:
                    best_model_metric = valid_corr
                    model.save_model(modelfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train APA-Net model")
    parser.add_argument(
        "--train_data", type=str, required=True, help="Path to training data file"
    )
    parser.add_argument(
        "--train_seq", type=str, required=True, help="Path to training sequences file"
    )
    parser.add_argument(
        "--valid_data", type=str, required=True, help="Path to validation data file"
    )
    parser.add_argument(
        "--valid_seq", type=str, required=True, help="Path to validation sequences file"
    )
    parser.add_argument(
        "--profiles", type=str, required=True, help="Path to cell type profiles file"
    )
    parser.add_argument(
        "--modelfile", type=str, required=True, help="Path to save the trained model"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=64,
        help="Batch size for training and validation",
    )
    parser.add_argument(
        "--epochs", type=int, default=200, help="Number of epochs for training"
    )
    parser.add_argument(
        "--project_name",
        type=str,
        default="APA-Net_Training",
        help="Wandb project name",
    )
    parser.add_argument(
        "--device", type=str, default="cuda:0", help="Device to run the training on"
    )
    parser.add_argument(
        "--use_wandb",
        type=str,
        default="True",
        help="Set to True to use wandb logging, False to disable it",
    )

    args = parser.parse_args()

    torch.manual_seed(7)
    random.seed(7)
    np.random.seed(7)

    train_data = np.load(args.train_data, allow_pickle=True)
    train_seq = np.load(args.train_seq, allow_pickle=True)
    valid_data = np.load(args.valid_data, allow_pickle=True)
    valid_seq = np.load(args.valid_seq, allow_pickle=True)
    profiles = pd.read_csv(args.profiles, index_col=0, sep="\t")

    config = {
        "batch_size": args.batch_size,
        "epochs": args.epochs,
        "project_name": args.project_name,
        "device": args.device,
        "opt": "Adam",
        "loss": "mse",
        "lr": 2.5e-05,
        "adam_weight_decay": 0.06,
        "conv1kc": 128,
        "conv1ks": 12,
        "conv1st": 1,
        "pool1ks": 25,
        "pool1st": 25,
        "cnvpdrop1": 0.2,
        "Matt_heads": 8,
        "Matt_drop": 0.2,
        "fc1_dims": [
            8192,
            4048,
            1024,
            512,
            256,
        ],  # first dimension will be calculated dynamically
        "fc1_dropouts": [0.3, 0.25, 0.25, 0.2, 0.1],
        "fc2_dims": [128, 32, 16, 1],  # first dimension will be calculated dynamically
        "fc2_dropouts": [0.2, 0.2, 0, 0],
    }

    main_train(
        train_seq,
        valid_seq,
        train_data,
        valid_data,
        profiles,
        args.modelfile,
        args.device,
        config,
        args.use_wandb,
    )
