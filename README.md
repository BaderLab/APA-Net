# APA-Net

APA-Net is a deep learning model designed for learning context specific APA usage. This guide covers the steps necessary to set up and run APA-Net.

## Installation

Before running APA-Net, ensure you have Python installed on your system. Clone this repository to your local machine:

```bash
git clone https://github.com/BaderLab/APA-Net.git
cd APA-Net

pip install .

```

# Usage

To train the APA-Net model, use the train_script.py script with the necessary command-line arguments:

```bash
python train_script.py \
--train_data "/path/to/train_data.npy" \
--train_seq "/path/to/train_seq.npy" \
--valid_data "/path/to/valid_data.npy" \
--valid_seq "/path/to/valid_seq.npy" \
--profiles "/path/to/celltype_profiles.tsv" \
--modelfile "/path/to/model_output.pt" \
--batch_size 64 \
--epochs 200 \
--project_name "APA-Net_Training" \
--device "cuda:1" \
--use_wandb "True"
```

# Arguments
- `--train_data`: Path to the training data file.
- `--train_seq`: Path to the training sequence data file.
- `--valid_data`: Path to the validation data file.
- `--valid_seq`: Path to the validation sequence data file.
- `--profiles`: Path to the cell type profiles file.
- `--modelfile`: Path where the trained model will be saved.
- `--batch_size`: Batch size for training (default: 64).
- `--epochs`: Number of training epochs (default: 200).
- `--project_name`: Name of the project for wandb logging.
- `--device`: Device to run the training on (e.g., 'cuda:1').
- `--use_wandb`: Flag to enable or disable wandb logging ('True' or 'False').

