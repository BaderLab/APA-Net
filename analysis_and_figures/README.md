# APA-Net Analysis and Figures

This directory contains all the analysis code, notebooks, and scripts used to generate the results and figures for the APA-Net research paper. The analysis pipeline covers the complete workflow from raw data processing to final figure generation.

## Quick Start

1. **Install Dependencies**:
   ```bash
   # Python packages
   pip install pandas numpy scipy matplotlib seaborn scikit-learn statsmodels jupyter
   
   # R packages
   Rscript -e "install.packages(c('dplyr', 'ggplot2', 'tidyr', 'viridis', 'patchwork', 'readxl', 'gridExtra', 'ggpubr', 'ggrepel', 'reshape2', 'corrplot', 'pheatmap', 'boot', 'Seurat', 'scCustomize'))"
   ```

2. **Run Core Analysis**:
   ```bash
   # Process data for APA-Net training
   cd data_processing && jupyter notebook Process_inputs_for_APA-Net.ipynb
   
   # Evaluate model performance  
   cd ../model_performance && jupyter notebook APA-NET_performance_plots.ipynb
   
   # Generate main figures
   cd ../visualization && jupyter notebook maaper_volcanos_barplots_figure6.ipynb
   ```

## Directory Structure

### 📊 `model_performance/`
Scripts for evaluating APA-Net model performance and interpreting learned features.

- **APA-NET_performance_plots.ipynb**: Model performance evaluation
  - Correlation analysis across cell types
  - Confidence intervals and statistical testing
  - Performance visualization plots
  
- **APA-Net_filter_interactions.ipynb**: Convolutional filter analysis
  - Extract filter activation patterns
  - Map filters to RNA-binding protein motifs
  - Analyze cell-type-specific patterns
  
- **APA-Net_heatmap_for_filter_interactions.ipynb**: Filter interaction heatmaps
  - Visualize RBP-filter relationships
  - Hierarchical clustering of binding patterns
  - Cell-type-specific interaction maps

### 🔧 `data_processing/`
Data preparation and preprocessing pipeline for APA-Net training.

- **Process_inputs_for_APA-Net.ipynb**: Main preprocessing pipeline
  - Sequence extraction and one-hot encoding
  - APA usage quantification and normalization
  - 5-fold cross-validation data splitting
  - Integration with cell-type-specific RBP profiles
  
- **APA_quantification_maaper_apalog_Dec2024.ipynb**: APA event quantification
  - MAAPER-based APA usage calculation
  - Statistical testing for differential usage
  - Quality control and filtering
  
- **emprical_fdr_thresholds_maaper_apalog.ipynb**: FDR threshold determination
  - Empirical null distribution generation
  - Cell-type-specific threshold calculation
  - Multiple testing correction validation

### 📈 `comparative_analysis/`
Comparative studies between APA changes and other molecular features.

- **APA_vs_DE.ipynb**: APA usage vs differential expression correlation
  - Spearman correlation analysis
  - Cell-type and condition-specific comparisons
  - Standardized effect size calculations
  
- **apa_correlation_across_celltypes.ipynb**: Cross-cell-type correlation analysis
  - APA usage pattern correlations between cell types
  - Shared and distinct regulatory mechanisms
  - Bootstrap confidence intervals
  
- **rbp_co_occurance_dissimilarity.ipynb**: RBP binding pattern analysis
  - Co-occurrence matrix construction
  - Dissimilarity analysis
  - Regulatory network inference

### 🎨 `visualization/`
Figure generation and visualization scripts.

- **maaper_volcanos_barplots_figure6.ipynb**: Main figure generation
  - Volcano plots for APA usage changes
  - Bar plots showing effect sizes by cell type
  - Statistical significance visualization
  - Multi-panel figure layout

### 🧬 `gene_expression/`
Differential gene expression analysis complementary to APA analysis.

- **DEG_ALS_genes.R**: ALS gene-specific expression analysis
- **DEG_MAST_analysis.R**: MAST-based differential expression
- **DEG_pathway_analysis.R**: Pathway enrichment for DEGs
- **DEG_visualization.R**: Expression data visualization

### 🛤️ `pathway_analysis/`
Gene set enrichment and pathway analysis for APA-affected genes.

- **APA_pathway_analysis.R**: Comprehensive pathway analysis
  - GO term enrichment analysis
  - Reactome pathway analysis
  - Custom gene set testing
  - Multiple testing correction

### 🔬 `preprocessing/`
Single-cell RNA-seq data preprocessing pipeline.

#### `processing_annotation/`
- **01_snRNA_cellranger_preprocess.sh**: Cell Ranger quantification
- **02_snRNA_process_QC.R**: Quality control and cell filtering
- **03_snRNA_clustering_annotation.R**: Cell clustering and type annotation
- **04a_snRNA_NSForest1.ipynb** & **04b_snRNA_NSForest2.ipynb**: NSForest cell type classification

#### `independent_datasets/`
- **01_read_matrices.R**: Independent dataset integration
- **02_harmony_int.R**: Batch correction with Harmony
- **03_doublet_removal_annotation.R**: Doublet detection and removal

### ⚙️ `functions/`
Utility functions and helper scripts.

- **export_gProfiler_results.R**: gProfiler result export utilities
- **readGMT.R**: GMT file reading functions
- **read_xlsx.R**: Excel file reading utilities
- **writeGMT.R**: GMT file writing functions

## Key Results Summary

### Model Performance
- **Overall correlation**: 0.585 (95% CI: 0.579-0.592)
- **Best cell type**: Microglia (r = 0.666)
- **Parameter count**: ~301M trainable parameters
- **Cross-validation**: 5-fold with sequence-based splitting

### Biological Insights
- **Condition correlations**: C9ALS vs sALS correlations range 0.65-0.84 across cell types
- **Cell-type specificity**: Distinct APA patterns per cell type
- **RBP associations**: 92 RBPs show significant filter associations
- **Pathway enrichment**: APA-affected genes enrich in neurodegeneration pathways

### Data Processing Stats
- **Samples processed**: ~57K training samples across 5 folds
- **Sequence length**: 4000 nucleotides (one-hot encoded)
- **Cell types**: 8 major brain cell types
- **Conditions**: Control, C9ALS, sALS

## Reproducibility

All analysis scripts are designed to be reproducible. Key considerations:

1. **Random seeds**: Set consistently across scripts (seed=1234)
2. **Package versions**: Use `renv` or `conda` for environment management
3. **Data paths**: Adjust file paths in notebooks to match your directory structure
4. **Computational resources**: Some analyses require substantial memory (>32GB recommended)

## Data Requirements

The analysis pipeline expects the following data files:
- Single-cell RNA-seq count matrices (10X format)
- APA usage quantification results (MAAPER output)
- Cell type annotations (CSV format)
- RBP expression profiles (TSV format)
- Reference genome and gene annotations

## Citation

If you use this analysis pipeline, please cite our work:
```bibtex
[Citation will be added upon publication]
```

## Support

For questions about specific analyses or to report issues:
1. Check script documentation and comments
2. Verify data file paths and formats
3. Ensure all dependencies are installed
4. Open an issue in the main repository

---

**Note**: Some analysis scripts may require adjustment of file paths to match your local data organization. The scripts assume a specific directory structure that may need to be adapted for your environment.