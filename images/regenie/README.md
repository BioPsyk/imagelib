# REGENIE Docker Image

REGENIE is a C++ program for whole genome regression modelling of large genome-wide association studies. It has been optimized to be computationally efficient for analyses of phenotypes with low trait prevalence.

**Note:** This image is only available for x86_64 (amd64) architecture as REGENIE binaries are not available for ARM64.

## Software Information

- **Software**: REGENIE
- **Version**: 3.6 (stable release)
- **License**: MIT
- **Repository**: https://github.com/rgcgithub/regenie
- **Documentation**: https://rgcgithub.github.io/regenie/

## Features

- Fast whole genome regression for quantitative and binary traits
- Handles population structure and relatedness
- Processes multiple phenotypes efficiently
- Supports gene/region-based tests and interaction tests
- Works with various genetic data formats (BGEN, PLINK, PLINK2)
- Memory-efficient implementation

## Usage

### Basic Usage

```bash
# Run with default help
docker run biopsyk/regenie:latest

# Run REGENIE step 1 (whole genome regression)
docker run -v /path/to/data:/data biopsyk/regenie:latest regenie \
  --step 1 \
  --bed /data/input \
  --phenoFile /data/phenotypes.txt \
  --covarFile /data/covariates.txt \
  --bsize 1000 \
  --out /data/output_step1

# Run REGENIE step 2 (association testing)
docker run -v /path/to/data:/data biopsyk/regenie:latest regenie \
  --step 2 \
  --bgen /data/input.bgen \
  --phenoFile /data/phenotypes.txt \
  --covarFile /data/covariates.txt \
  --pred /data/output_step1_pred.list \
  --bsize 400 \
  --out /data/output_step2
```

### Common Options

- `--step`: Specify step 1 (whole genome regression) or step 2 (association testing)
- `--bed/--bgen/--pgen`: Input genetic data files
- `--phenoFile`: Phenotype file
- `--covarFile`: Covariate file
- `--bsize`: Block size for computation
- `--out`: Output file prefix
- `--threads`: Number of threads to use

### Example with Multiple Phenotypes

```bash
docker run -v /path/to/data:/data biopsyk/regenie:latest regenie \
  --step 1 \
  --bed /data/input \
  --phenoColList pheno1,pheno2,pheno3 \
  --phenoFile /data/phenotypes.txt \
  --covarFile /data/covariates.txt \
  --bsize 1000 \
  --lowmem \
  --out /data/multi_pheno_step1
```

### Binary Traits Analysis

```bash
docker run -v /path/to/data:/data biopsyk/regenie:latest regenie \
  --step 2 \
  --bgen /data/input.bgen \
  --phenoFile /data/phenotypes.txt \
  --covarFile /data/covariates.txt \
  --bt \
  --firth --approx \
  --pred /data/step1_pred.list \
  --bsize 400 \
  --out /data/binary_trait_results
```

## Environment Variables

No specific environment variables are required.

## Data Format Requirements

### Phenotype File
Tab or space-delimited file with:
- Column 1: Family ID
- Column 2: Individual ID
- Column 3+: Phenotype values

### Covariate File
Tab or space-delimited file with:
- Column 1: Family ID
- Column 2: Individual ID
- Column 3+: Covariate values

## Performance Tips

1. Use `--lowmem` option for large datasets to reduce memory usage
2. Adjust `--bsize` based on available memory (larger blocks = faster but more memory)
3. Use `--threads` to parallelize computation
4. For binary traits with low case counts, use `--firth` for better convergence

## Citation

If you use REGENIE in your research, please cite:

Mbatchou, J., Barnard, L., Backman, J. et al. Computationally efficient whole-genome regression for quantitative and binary traits. Nat Genet 53, 1097–1103 (2021). https://doi.org/10.1038/s41588-021-00870-7

## Support

For issues or questions about REGENIE:
- GitHub Issues: https://github.com/rgcgithub/regenie/issues
- Documentation: https://rgcgithub.github.io/regenie/

For Docker image specific issues, please contact your institute's bioinformatics support.