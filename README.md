# TMScal: Tumor Mutational Signature Calculator

## Overview
TMScal is an R package for predicting prognostic risks using Tumor Mutational 
Spectrum (TMS) from whole exome sequencing (WES) data of MSS (microsatellite 
stable) cancer patients.

The method uses 192 trinucleotide mutational pathways with adaptive 
stability-based feature selection.

## Key Features
- Data preparation (TCGA format conversion, MSS filtering)
- Pathway stability analysis (bootstrap-based)
- Feature selection (stability cutoff 0.07 + adaptive P-value)
- Model training (60:40 stratified split)
- Independent validation
- Prognostic risk prediction (Low/Medium/High)
- Visualization (KM curves, stability plots, feature importance)

## Installation
```r
# Install from GitHub
devtools::install_github("cxue/TMScal")
```

## Quick Start
```r
library(TMScal)

# One-click pipeline
result <- run_all_pipeline(
    maf_file = "data_mutations.txt",
    clin_file = "data_clinical.txt",
    output_dir = "./tmscal_output"
)

# View results
summary(result$model)
plot_km_curve(result$model, result$clin_data)
```

## Data Requirements
- WES (whole exome sequencing) data
- MSS (microsatellite stable) samples
- Minimum: 300 samples, 50 events
- Recommended: 400+ samples, 70+ events

## Methodology
- Extract 192 trinucleotide mutational pathways
- CLR transformation of pathway proportions
- Bootstrap-based stability assessment (200 iterations)
- Univariate Cox regression for prognostic value
- Feature selection: stability < 0.07 + adaptive P-value cutoff
- 60:40 stratified train/test split
- TMS score = weighted sum of selected pathway CLR values
- Tertile-based risk stratification (Low/Medium/High)

## LICENSE：

  Apache License (>= 2.0) + file LICENSE
