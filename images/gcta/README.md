# GCTA Docker Image

This Docker image provides GCTA (Genome-wide Complex Trait Analysis), a software package for the analysis of data from genome-wide association studies (GWASs).

## Description

GCTA (Genome-wide Complex Trait Analysis) was originally designed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for complex traits (the GREML method), and has subsequently been extended for many other analyses to better understand the genetic architecture of complex traits.

## Version Information

- Internal image version: 1.0.0
- GCTA version: 1.94.1 (can be verified by running `gcta64 --help`)

## Key Features

- **GREML**: Estimating SNP-based heritability
- **COJO**: Conditional & joint association analysis using GWAS summary statistics
- **fastGWA**: Ultra-fast genome-wide association analysis
- **MLMA**: Mixed linear model association analysis
- **GRM**: Genetic relationship matrix estimation
- **PCA**: Principal component analysis
- **GSMR**: Generalized summary-data-based Mendelian randomization

## Usage

### Basic Usage

```bash
# Estimate genetic relationship matrix (GRM)
docker run -v $(pwd):/data biopsyk/gcta gcta64 --bfile mydata --autosome --make-grm --out test

# GREML analysis to estimate heritability
docker run -v $(pwd):/data biopsyk/gcta gcta64 --reml --grm test --pheno mydata.phen --out test

# FastGWA analysis
docker run -v $(pwd):/data biopsyk/gcta gcta64 --fastGWA --bfile mydata --pheno mydata.phen --out fastgwa_result

# COJO analysis using summary statistics
docker run -v $(pwd):/data biopsyk/gcta gcta64 --cojo-file gwas_summary.txt --ref-ld-chr ld_ref --out cojo_result
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/gcta:1.94.1

# Run GREML analysis
singularity exec gcta_1.94.1.sif gcta64 --reml --grm test --pheno mydata.phen --out test
```

## Input Formats

### Genotype Data
- PLINK binary format (.bed, .bim, .fam)
- PLINK text format (.ped, .map)
- Dosage format (MACH/minimac output)

### Phenotype Data
- Plain text files with family ID, individual ID, and phenotype values
- Missing values should be coded as -9 or NA

### Summary Statistics
- Tab-delimited text files with SNP, chromosome, position, effect allele, other allele, frequency, beta, SE, P-value

## Common Analysis Workflows

### 1. SNP-based Heritability Estimation

```bash
# Step 1: Create GRM
docker run -v $(pwd):/data biopsyk/gcta gcta64 --bfile mydata --autosome --maf 0.01 --make-grm --out mydata_grm

# Step 2: GREML analysis
docker run -v $(pwd):/data biopsyk/gcta gcta64 --reml --grm mydata_grm --pheno mydata.phen --out heritability_result
```

### 2. Genome-wide Association Analysis

```bash
# FastGWA (recommended for large datasets)
docker run -v $(pwd):/data biopsyk/gcta gcta64 --fastGWA --bfile mydata --pheno mydata.phen --thread-num 4 --out gwas_result

# Traditional MLMA
docker run -v $(pwd):/data biopsyk/gcta gcta64 --mlma --bfile mydata --grm mydata_grm --pheno mydata.phen --out mlma_result
```

### 3. Conditional Analysis with COJO

```bash
# Identify independent signals
docker run -v $(pwd):/data biopsyk/gcta gcta64 --cojo-file gwas_summary.ma --ref-ld-chr reference_data --cojo-slct --out cojo_signals
```

## Performance Notes

- Use `--thread-num` to specify number of CPU cores
- For large datasets (>50K individuals), consider using fastGWA instead of MLMA
- GRM computation is memory-intensive; use `--make-grm-part` for very large samples

## Examples

```bash
# Multi-step heritability analysis with covariates
docker run -v $(pwd):/data biopsyk/gcta gcta64 --bfile mydata --autosome --maf 0.01 --make-grm --out mydata_grm
docker run -v $(pwd):/data biopsyk/gcta gcta64 --reml --grm mydata_grm --pheno mydata.phen --qcovar pcs.txt --out h2_with_pcs

# Estimate genetic correlation between two traits
docker run -v $(pwd):/data biopsyk/gcta gcta64 --reml-bivar --grm mydata_grm --pheno two_traits.phen --out genetic_correlation
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- GCTA requires substantial memory for large datasets
- Use absolute paths or ensure all input files are in the mounted directory

## Memory Requirements

| Analysis Type | Sample Size | Approximate Memory |
|---------------|-------------|-------------------|
| GRM creation | 10K | ~4 GB |
| GRM creation | 50K | ~100 GB |
| GREML | 10K | ~8 GB |
| FastGWA | 100K+ | ~16 GB |

## Building the Image

From this directory:

```bash
docker build -t biopsyk/gcta:latest .
```

## References

- [GCTA Software Website](https://yanglab.westlake.edu.cn/software/gcta/)
- [GCTA GitHub Repository](https://github.com/jianyangqt/gcta)
- [GCTA Documentation](https://yanglab.westlake.edu.cn/software/gcta/index.html)
- Yang et al. (2011) GCTA: a tool for Genome-wide Complex Trait Analysis. Am J Hum Genet. 88(1): 76-82. 