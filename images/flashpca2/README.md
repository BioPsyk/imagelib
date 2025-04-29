# FlashPCA2 Docker Image

This Docker image provides FlashPCA2 (version 2.1), a fast implementation of principal component analysis (PCA) for large-scale genome-wide data.

## Description

FlashPCA2 is a fast implementation of principal component analysis (PCA) for large-scale genome-wide data. It is particularly useful for analyzing genetic data and can handle large datasets efficiently.

## Version Information

- Internal image version: 1.0.0
- FlashPCA2 version: 2.1 (can be verified by running `flashpca --version`)

## Usage

### Basic Usage

```bash
# Run flashpca2 on a PLINK binary file
docker run -v $(pwd):/data biopsyk/flashpca2 --bfile your_data

# Run with specific number of components
docker run -v $(pwd):/data biopsyk/flashpca2 --bfile your_data --ndim 10
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/flashpca2:2.1

# Run flashpca2
singularity exec flashpca2_2.1.sif flashpca --bfile your_data
```

## Input/Output Files

### Input
- PLINK binary files (.bed, .bim, .fam)
- Other supported formats as per flashpca2 documentation

### Output
By default, flashpca2 produces the following files:
- `eigenvectors.txt`: Top k eigenvectors
- `pcs.txt`: Top k principal components
- `eigenvalues.txt`: Top k eigenvalues
- `pve.txt`: Proportion of variance explained

## Example

```bash
# Assuming you have PLINK binary files in the current directory
docker run -v $(pwd):/data biopsyk/flashpca2 --bfile your_data --ndim 10

# This will produce:
# - eigenvectors.txt
# - pcs.txt
# - eigenvalues.txt
# - pve.txt
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- The container includes only the necessary runtime dependencies to minimize size

## References

- [FlashPCA2 GitHub Repository](https://github.com/gabraham/flashpca)
- [FlashPCA2 Documentation](https://github.com/gabraham/flashpca/blob/master/README.md) 