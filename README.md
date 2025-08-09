# APA-Net

APA-Net is a deep learning model designed for learning context-specific APA (Alternative Polyadenylation) usage. This guide covers the steps necessary to set up and run APA-Net.

## Requirements

- Python 3.8 or higher
- PyTorch 1.8.0 or higher
- NumPy
- Pandas
- SciPy
- tqdm
- wandb (optional, for experiment tracking)

## Installation

### Option 1: Install from source (Recommended)

1. Clone this repository to your local machine:
```bash
git clone https://github.com/BaderLab/APA-Net.git
cd APA-Net
```

2. Install dependencies manually for better control:
```bash
# For CPU-only version (smaller download)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# For GPU version (if you have CUDA)
pip install torch torchvision torchaudio

# Install other dependencies
pip install numpy pandas scipy tqdm wandb
```

3. Install the package:
```bash
pip install .
```

### Option 2: One-command installation
```bash
pip install .
```
*Note: This will install the full PyTorch with CUDA support, which is a large download (~2GB).*

## Data Format

APA-Net expects input data in `.npy` format with the following structure:
- **Shape**: `(n_samples, 9)` where each row represents one sample
- **Columns**:
  - Column 0: Float value (sample ID/index)
  - Column 1: String (cell type name)
  - Column 2: String (additional metadata)
  - Column 3: Float value 
  - Column 4: String (additional metadata)
  - Column 5: String (genomic coordinates/switch name)
  - Column 6: NumPy array of shape `(4, 4000)` - one-hot encoded DNA sequence
  - Column 7: Float (target APA usage value)
  - Column 8: NumPy array of shape `(327,)` - cell type profile features

## Usage

### Training the Model

To train the APA-Net model, use the train_script.py script:

```bash
cd apamodel
python train_script.py \
  --train_data "/path/to/train_data.npy" \
  --valid_data "/path/to/valid_data.npy" \
  --modelfile "/path/to/model_output.pt" \
  --batch_size 64 \
  --epochs 200 \
  --device "cpu" \
  --use_wandb "False"
```

### Testing the Model

You can test the model with sample data:

```bash
# Create a simple test script
python -c "
import sys
sys.path.append('./apamodel')
from model import APANET, APAData
import numpy as np
import torch

# Load your data
data = np.load('your_data.npy', allow_pickle=True)

# Configure model (using CPU)
config = {
    'device': 'cpu',
    'opt': 'Adam',
    'loss': 'mse',
    'lr': 2.5e-05,
    'adam_weight_decay': 0.09,
    'conv1kc': 128,
    'conv1ks': 12,
    'conv1st': 1,
    'pool1ks': 16,
    'pool1st': 16,
    'cnvpdrop1': 0,
    'Matt_heads': 8,
    'Matt_drop': 0.2,
    'fc1_dims': [8192, 4048, 1024, 512, 256],
    'fc1_dropouts': [0.25, 0.25, 0.25, 0, 0],
    'fc2_dims': [128, 32, 16, 1],
    'fc2_dropouts': [0.2, 0.2, 0, 0],
    'psa_query_dim': 128,
    'psa_num_layers': 1,
    'psa_nhead': 1,
    'psa_dim_feedforward': 1024,
    'psa_dropout': 0
}

# Create and test model
model = APANET(config)
model.compile()
print('Model created successfully!')
"
```

## Command Line Arguments

- `--train_data`: Path to the training data file (required)
- `--valid_data`: Path to the validation data file (required)
- `--modelfile`: Path where the trained model will be saved (required)
- `--batch_size`: Batch size for training (default: 64)
- `--epochs`: Number of training epochs (default: 200)
- `--project_name`: Name of the project for wandb logging (default: "APA-Net_Training")
- `--device`: Device to run the training on - use "cpu" or "cuda:0" (default: "cuda:0")
- `--use_wandb`: Enable wandb logging - "True" or "False" (default: "True")

## Model Architecture

APA-Net is a deep neural network that combines:
- **Convolutional layers** for sequence feature extraction
- **Self-attention mechanism** for capturing long-range dependencies
- **Fully connected layers** for prediction
- **Cell type profile integration** for context-specific modeling

The model has approximately 301M parameters and processes:
- Input: DNA sequences (4×4000) + cell type profiles (327 features)
- Output: APA usage prediction (single value)

## Troubleshooting

### Common Issues

1. **CUDA errors**: If you encounter CUDA-related errors, install the CPU-only version of PyTorch:
   ```bash
   pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
   ```

2. **Memory issues**: Reduce batch size if you encounter out-of-memory errors:
   ```bash
   --batch_size 32
   ```

3. **Data format errors**: Ensure your data has the correct shape `(n_samples, 9)` with sequences of shape `(4, 4000)` and cell type profiles of shape `(327,)`.

### CPU vs GPU Usage

- **CPU**: Slower but more compatible. Use `--device "cpu"`
- **GPU**: Faster training. Use `--device "cuda:0"` (requires CUDA-compatible PyTorch installation)

## Example

Here's a complete example of training APA-Net:

```bash
# Navigate to the model directory
cd APA-Net/apamodel

# Train the model
python train_script.py \
  --train_data "../test_fold_0.npy" \
  --valid_data "../test_fold_0.npy" \
  --modelfile "./trained_model.pt" \
  --batch_size 32 \
  --epochs 50 \
  --device "cpu" \
  --use_wandb "False" \
  --project_name "APA-Net_Test"
```

## Analysis and Figures

The `analysis_and_figures/` directory contains all the code and notebooks used to reproduce the results and figures from our APA-Net research paper. This comprehensive analysis pipeline covers data processing, model evaluation, comparative analysis, and visualization.

### Directory Structure

```
analysis_and_figures/
├── model_performance/          # APA-Net model evaluation and performance analysis
├── data_processing/            # Data preparation and preprocessing for APA-Net
├── comparative_analysis/       # Comparative studies (APA vs DE, correlations)
├── visualization/              # Figure generation and plotting scripts
├── gene_expression/            # Differential gene expression analysis
├── pathway_analysis/           # Gene set enrichment and pathway analysis
├── preprocessing/              # Single-cell RNA-seq data preprocessing pipeline
└── functions/                  # Utility functions and helper scripts
```

### Getting Started with Analysis

1. **Prerequisites**: Make sure you have the following R and Python packages installed:

**R packages:**
```r
install.packages(c("dplyr", "ggplot2", "tidyr", "viridis", "patchwork", 
                   "readxl", "gridExtra", "ggpubr", "ggrepel", "reshape2", 
                   "corrplot", "pheatmap", "boot", "Seurat", "scCustomize"))
```

**Python packages:**
```bash
pip install pandas numpy scipy matplotlib seaborn scikit-learn statsmodels
```

2. **Data Requirements**: The analysis scripts expect data in specific locations. You may need to adjust file paths in the notebooks to match your data directory structure.

### Analysis Modules

#### 1. Model Performance (`model_performance/`)
- **APA-NET_performance_plots.ipynb**: Generates correlation plots showing model performance across cell types
- **APA-Net_filter_interactions.ipynb**: Analyzes convolutional filter interactions and RBP binding patterns  
- **APA-Net_heatmap_for_filter_interactions.ipynb**: Creates heatmaps showing filter-RBP interactions

#### 2. Data Processing (`data_processing/`)
- **Process_inputs_for_APA-Net.ipynb**: Main data preprocessing pipeline for APA-Net training data
  - Processes RNA sequences and APA usage data
  - Generates one-hot encoded sequences
  - Creates 5-fold cross-validation splits
  - Formats data for model training
- **APA_quantification_maaper_apalog_Dec2024.ipynb**: APA event quantification using MAAPER
- **emprical_fdr_thresholds_maaper_apalog.ipynb**: Determines empirical FDR thresholds for significance testing

#### 3. Comparative Analysis (`comparative_analysis/`)
- **APA_vs_DE.ipynb**: Compares APA changes with differential expression
  - Correlation analysis between APA usage and gene expression changes
  - Cell-type-specific comparisons
  - Statistical significance testing
- **apa_correlation_across_celltypes.ipynb**: Cross-cell-type APA correlation analysis
- **rbp_co_occurance_dissimilarity.ipynb**: RNA-binding protein co-occurrence analysis

#### 4. Visualization (`visualization/`)
- **maaper_volcanos_barplots_figure6.ipynb**: Creates volcano plots and bar plots for Figure 6
  - APA usage changes across conditions
  - Cell-type-specific visualizations
  - Statistical significance visualization

#### 5. Gene Expression (`gene_expression/`)
- **DEG_ALS_genes.R**: Analysis of ALS-associated gene expression
- **DEG_MAST_analysis.R**: MAST-based differential expression analysis
- **DEG_pathway_analysis.R**: Pathway enrichment analysis for DEGs
- **DEG_visualization.R**: Visualization of differential expression results

#### 6. Pathway Analysis (`pathway_analysis/`)
- **APA_pathway_analysis.R**: Gene set enrichment analysis for APA-affected genes
  - GO term enrichment
  - Reactome pathway analysis
  - Custom gene set analysis

#### 7. Preprocessing (`preprocessing/`)
- **processing_annotation/**: Single-cell RNA-seq processing pipeline
  - `01_snRNA_cellranger_preprocess.sh`: Cell Ranger preprocessing
  - `02_snRNA_process_QC.R`: Quality control and filtering
  - `03_snRNA_clustering_annotation.R`: Cell clustering and annotation
  - `04a_snRNA_NSForest1.ipynb` & `04b_snRNA_NSForest2.ipynb`: NSForest cell type classification
- **independent_datasets/**: Processing of additional validation datasets
  - `01_read_matrices.R`: Matrix reading and preprocessing
  - `02_harmony_int.R`: Harmony integration for batch correction
  - `03_doublet_removal_annotation.R`: Doublet detection and removal

### Reproducing Key Results

#### Figure Generation
To reproduce the main figures from the paper:

1. **Model Performance Plots**:
   ```bash
   cd analysis_and_figures/model_performance
   jupyter notebook APA-NET_performance_plots.ipynb
   ```

2. **APA Usage Analysis**:
   ```bash
   cd analysis_and_figures/visualization  
   jupyter notebook maaper_volcanos_barplots_figure6.ipynb
   ```

3. **Comparative Analysis**:
   ```bash
   cd analysis_and_figures/comparative_analysis
   jupyter notebook APA_vs_DE.ipynb
   ```

#### Data Processing Pipeline
To process your own data through the complete pipeline:

1. **Start with raw single-cell data**:
   ```bash
   cd analysis_and_figures/preprocessing/processing_annotation
   bash 01_snRNA_cellranger_preprocess.sh
   ```

2. **Process and prepare for APA-Net**:
   ```bash
   cd analysis_and_figures/data_processing
   jupyter notebook Process_inputs_for_APA-Net.ipynb
   ```

### Key Results and Interpretations

- **Model Performance**: APA-Net achieves correlation coefficients of 0.56-0.67 across cell types
- **Cell-Type Specificity**: Microglia show highest model performance, indicating stronger APA regulatory patterns
- **Condition Comparison**: Strong correlations (0.65-0.84) between C9ALS and sALS APA changes across cell types
- **Biological Validation**: APA changes correlate with known ALS pathways and RBP targets

### Data Availability

The analysis scripts reference several data sources:
- Single-cell RNA-seq count matrices
- APA usage quantification results  
- Cell type annotations
- RBP expression profiles
- Reference genome and annotations

Please ensure you have access to the appropriate datasets before running the analysis scripts.

### Citation

If you use this analysis pipeline, please cite our paper:
```
[Paper citation to be added upon publication]
```

For questions about the analysis pipeline, please open an issue in the GitHub repository.

