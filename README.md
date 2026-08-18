# TMScal: Tumor Mutational Signature Calculator

## Overview
TMScal is an R package for predicting prognostic risks using Tumor Mutational 
Spectrum (TMS) from whole exome sequencing (WES) data of MSS (microsatellite 
stable) cancer patients.

The method uses 192 trinucleotide mutational pathways with adaptive 
feature selection.

## Key Features
- Data preparation (TCGA format conversion, MSS filtering)
- Model training with internal CV tuning (frequency threshold, alpha, lambda)
- Prognostic risk prediction (Low/Medium/High)
- Visualization (KM curves)

## Installation
```r
# Install from GitHub
devtools::install_github("cxue/TMScal")
```

## Quick Start
```r
library(TMScal)

# Example: SYSUCC COADREAD MSS data (included)
maf_file <- system.file("extdata", "data_mutations_train_70.txt", package = "TMScal")
clin_file <- system.file("extdata", "data_clinical_train_70.txt", package = "TMScal")

# Train model
model <- train(maf_file, clin_file, output_dir = "./model")

# View results
summary(model)

# Predict test set
test_maf <- system.file("extdata", "data_mutations_test_30.txt", package = "TMScal")
predictions <- predict(test_maf, model, output_file = "./predictions.csv")

# Plot KM curve
plot_km_curve("./predictions.csv", 
              system.file("extdata", "data_clinical_test_30.txt", package = "TMScal"),
              output_file = "./KM_curve.png")
			  
```

## Data Requirements
- WES (whole exome sequencing) data
- MSS (microsatellite stable) samples
- Minimum: 300 samples, 50 events
- Recommended: 400+ samples, 70+ events
- MAF columns: chromosome, start_pos, ref_allele, alt_allele, sample_id, variant_class, context (see extdata/data_mutations_train_70.txt as an example)
- Clinical file columns: sample_id, os_time, os_event (see extdata/data_clinical_train_70.txt as an example)

## Data Format
MAF file columns:
```text
chromosome  start_pos  ref_allele  alt_allele  sample_id  variant_class  context
```

Clinical file columns:
```text
sample_id  os_time  os_event
```

## Included Example Data
The package includes SYSUCC COADREAD MSS data:
- data_mutations_train_70.txt - Training MAF (483 samples)
- data_clinical_train_70.txt - Training clinical data
- data_mutations_test_30.txt - Test MAF (206 samples)
- data_clinical_test_30.txt - Test clinical data
- data_mutations.SYSUCC-COADRECT_MSS.matched.txt - Full MAF (689 samples)
- data_clinical.SYSUCC-COADRECT_MSS.matched.txt - Full clinical data

## Methodology
- Extract 192 trinucleotide mutational pathways from MAF
- Build pathway proportion matrix (Pi = count / total mutations)
- Internal 5-fold CV for frequency threshold selection (0.05-0.15)
- Internal 5-fold CV for alpha selection (0.05-0.50)
- Lambda search (10-60 features range)
- Elastic Net Cox regression (log10(Pi+1) + TMB)
- 70:30 stratified train/test split
- TMS score = weighted sum of selected pathway values
- Tertile-based risk stratification (Low/Medium/High)

Performance (SYSUCC COADREAD MSS)
- Training: 483 samples, 99 events
- Test: 206 samples, 42 events
- Test HR (High vs Low): 2.18 (1.04-4.58)
- Test HR P: 0.04
- Test C-index: 0.579


## LICENSE：

  Apache License (>= 2.0) + file LICENSE
  The methods and systems described herein are the subject of a pending
PCT international patent application (Application No.: PCT/CN2026/112587,
filed on July 22, 2026).