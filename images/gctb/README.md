# GCTB Docker Image

This Docker image provides GCTB (Genome-wide Complex Trait Bayesian analysis), a software tool that comprises a family of Bayesian linear mixed models for complex trait analyses using genome-wide SNPs.

## Description

GCTB was developed to simultaneously estimate the joint effects of all SNPs and the genetic architecture parameters for a complex trait, including SNP-based heritability, polygenicity and the joint distribution of effect sizes and minor allele frequencies. Version 2.0+ includes summary-data-based versions of the individual-level data Bayesian linear mixed models.

## Version Information

- Internal image version: 1.0.0
- GCTB version: 2.01 Beta (can be verified by running `gctb`)

## Key Features

- **SBayesR**: Summary-data-based Bayesian regression for polygenic prediction
- **SBayesS**: Summary-data-based sparse Bayesian regression
- **SBayesC**: Summary-data-based Bayesian regression with continuous shrinkage
- **SBayesRC**: Multi-component summary-data-based Bayesian regression
- **BayesN**: Individual-level data Bayesian regression models
- **Cross-validation**: K-fold cross-validation for model evaluation
- **GWAS**: Genome-wide association analysis using Bayesian methods

## Usage

### Basic Usage

```bash
# SBayesR analysis using GWAS summary statistics
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_ref.ldm.sparse \
    --pi 0.95,0.02,0.02,0.01 \
    --gamma 0.0,0.01,0.1,1 \
    --gwas-summary gwas_summary.ma \
    --chain-length 10000 \
    --burn-in 2000 \
    --out-freq 10 \
    --out sbayesr_result

# SBayesS sparse regression analysis
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes S \
    --ldm ld_ref.ldm.sparse \
    --gwas-summary gwas_summary.ma \
    --chain-length 10000 \
    --hsq 0.5 \
    --out sbayess_result

# Cross-validation analysis
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_ref.ldm.sparse \
    --gwas-summary gwas_summary.ma \
    --cross-validation \
    --num-folds 5 \
    --out cv_result
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/gctb:2.05

# Run SBayesR analysis
singularity exec gctb_2.05.sif gctb --sbayes R --mldm ld_ref.ldm.sparse --gwas-summary gwas_summary.ma --out result
```

## Input Formats

### GWAS Summary Statistics
- Tab-delimited text files with columns: SNP, A1, A2, freq, b, se, p, N
- Compressed files (.gz) are supported

### LD Reference Data
- Sparse LD matrix format (.ldm.sparse)
- Dense LD matrix format (.ldm.dense)
- GCTA binary LD format

### Individual-level Data
- PLINK binary format (.bed, .bim, .fam)
- PLINK text format (.ped, .map)

## Common Analysis Workflows

### 1. SBayesR for Polygenic Risk Scoring

```bash
# Step 1: Prepare LD reference (if needed)
docker run -v $(pwd):/data biopsyk/gctb gctb --bfile reference_panel \
    --make-sparse-ldm \
    --chisq-max 80 \
    --out ld_reference

# Step 2: Run SBayesR
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_reference.ldm.sparse \
    --pi 0.95,0.02,0.02,0.01 \
    --gamma 0.0,0.01,0.1,1 \
    --gwas-summary gwas_summary.ma \
    --chain-length 10000 \
    --burn-in 2000 \
    --out sbayesr_prs
```

### 2. Model Comparison with Cross-validation

```bash
# Compare SBayesR vs SBayesS
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_ref.ldm.sparse \
    --gwas-summary gwas_summary.ma \
    --cross-validation \
    --num-folds 5 \
    --out sbayesr_cv

docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes S \
    --ldm ld_ref.ldm.sparse \
    --gwas-summary gwas_summary.ma \
    --cross-validation \
    --num-folds 5 \
    --out sbayess_cv
```

### 3. Multi-component Analysis

```bash
# SBayesRC with custom priors
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes RC \
    --mldm ld_ref.ldm.sparse \
    --annot functional_annotations.txt \
    --gwas-summary gwas_summary.ma \
    --chain-length 15000 \
    --burn-in 3000 \
    --out sbayesrc_result
```

## Model Parameters

### SBayesR Parameters
- `--pi`: Mixture proportions for effect size categories
- `--gamma`: Variance scaling factors for effect size categories
- `--hsq`: SNP-based heritability (if not estimated)

### SBayesS Parameters
- `--pi`: Proportion of non-zero effects
- `--lambda`: Regularization parameter
- `--hsq`: SNP-based heritability

### MCMC Parameters
- `--chain-length`: Number of MCMC iterations
- `--burn-in`: Number of burn-in iterations
- `--out-freq`: Output frequency for MCMC samples

## Output Files

### Main Results
- `.parRes`: Parameter estimates and model fit statistics
- `.snpRes`: SNP-wise posterior effect sizes and inclusion probabilities
- `.mcmcsamples`: MCMC samples of model parameters

### Cross-validation Results
- `.cv`: Cross-validation prediction accuracies
- `.pred`: Individual-level predictions (when applicable)

## Performance Notes

- Use sparse LD matrices for faster computation
- Increase `--chain-length` for better convergence
- Monitor MCMC traces to ensure proper mixing
- Use multiple chains for convergence diagnostics

## Examples

```bash
# Quick SBayesR analysis with default parameters
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_ref.ldm.sparse \
    --gwas-summary gwas_summary.ma \
    --out quick_sbayesr

# Advanced analysis with custom priors and longer chain
docker run -v $(pwd):/data biopsyk/gctb gctb --sbayes R \
    --mldm ld_ref.ldm.sparse \
    --pi 0.99,0.005,0.004,0.001 \
    --gamma 0.0,0.001,0.01,0.1 \
    --gwas-summary gwas_summary.ma \
    --chain-length 20000 \
    --burn-in 5000 \
    --out detailed_sbayesr
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- GCTB analyses can be computationally intensive
- Ensure sufficient memory for large datasets

## Memory Requirements

| Analysis Type | SNPs | Approximate Memory |
|---------------|------|-------------------|
| SBayesR | 1M | ~8 GB |
| SBayesR | 5M | ~32 GB |
| Cross-validation | 1M | ~16 GB |

## Building the Image

From this directory:

```bash
docker build -t biopsyk/gctb:latest .
```

## References

- [GCTB Software Website](https://cnsgenomics.com/software/gctb/)
- [GCTB GitHub Repository](https://github.com/jianyangqt/gctb)
- Lloyd-Jones et al. (2019) Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nature Communications. 10:5086
- Zeng et al. (2021) Signatures of negative selection in the genetic architecture of human complex traits. Nature Genetics. 53:746-754 