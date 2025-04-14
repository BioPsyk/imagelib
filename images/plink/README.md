# PLINK Docker Image

This Docker image contains both PLINK 1.9 and PLINK 2, commonly used tools for whole genome association analysis.

## Included Software

- PLINK 1.9 (latest stable)
- PLINK 2 (alpha 4.0)

## Usage

The image is designed to be used with data mounted from the host system. The working directory is set to `/data`.

### Basic Usage

```bash
# Run PLINK 1.9
docker run -v $(pwd):/data biopsyk/plink:latest plink --help

# Run PLINK 2
docker run -v $(pwd):/data biopsyk/plink:latest plink2 --help

# Run analysis with mounted data
docker run -v $(pwd):/data biopsyk/plink:latest plink --bfile mydata --assoc --out results
```

### Example Commands

1. Basic association analysis:
```bash
docker run -v $(pwd):/data biopsyk/plink:latest plink --bfile mydata --assoc --out results
```

2. Quality control:
```bash
docker run -v $(pwd):/data biopsyk/plink:latest plink2 \
    --bfile mydata \
    --maf 0.01 \
    --geno 0.1 \
    --hwe 1e-6 \
    --make-bed \
    --out qc_data
```

## Building the Image

From this directory:

```bash
docker build -t biopsyk/plink:latest .
```

## Version Information

The versions of installed software can be found in `/etc/plink-version` inside the container. 