#!/bin/bash

# Define common variables
TRANSCRIPTOME="/ALS_snRNA/refdata-gex-GRCh38-2020-A"
FASTQ_DIR="/ALS_snRNA/fastq/"
EXPECTED_CELLS=6000

# Function to run cellranger count
run_cellranger() {
    local sample_id=$1
    
    cellranger count --id="${sample_id}_snRNA" \
                     --transcriptome="$TRANSCRIPTOME" \
                     --fastqs="$FASTQ_DIR" \
                     --include-introns \
                     --sample="${sample_id}_snRNA" \
                     --expect-cells="$EXPECTED_CELLS"
}

# Arrays of sample prefixes
CTRL_SAMPLES=(CTRL{1..6})
C9ALS_NO_FTLD_SAMPLES=(C9ALSnoFTLD{1..3})
C9ALS_FTLD_SAMPLES=(C9ALSFTLD{1..6})
SALS_NO_FTLD_SAMPLES=(sALSnoFTLD{1..8})

# Process CTRL samples
for sample in "${CTRL_SAMPLES[@]}"; do
    run_cellranger "$sample"
done

# Process C9ALS without FTLD samples
for sample in "${C9ALS_NO_FTLD_SAMPLES[@]}"; do
    run_cellranger "$sample"
done

# Process C9ALS with FTLD samples
for sample in "${C9ALS_FTLD_SAMPLES[@]}"; do
    run_cellranger "$sample"
done

# Process sALS without FTLD samples
for sample in "${SALS_NO_FTLD_SAMPLES[@]}"; do
    run_cellranger "$sample"
done