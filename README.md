# Institute Docker Images

This repository contains Dockerfiles for commonly used software at our institute. Each image is optimized for size and built for multiple architectures (amd64 and arm64).

## Available Images

| Image Name | Description | Size | Latest Tag | Architectures | Software Version | Image Version |
|------------|-------------|------|------------|---------------|------------------|---------------|
| plink | PLINK1.9 and PLINK2 for genetic analysis | ~100MB | latest | amd64, arm64 | 1.9 & 2.0 | 1.0.0 |
| flashpca2 | Fast PCA implementation for genome-wide data | ~100MB | latest | amd64, arm64 | 2.1 | 1.0.0 |
| gcta | Genome-wide Complex Trait Analysis | ~150MB | latest | amd64, arm64 | 1.94.1 | 1.0.0 |
| gctb | Bayesian linear mixed models for complex traits | ~120MB | latest | amd64, arm64 | 2.05 | 1.0.0 |
| readstat | Statistical data format conversion tool | ~120MB | latest | amd64, arm64 | 1.1.9 | 1.0.0 |
| regenie | Whole genome regression for GWAS | ~150MB | latest | amd64 | 3.6 | 1.0.0 |
| bcftools | VCF/BCF manipulation with HTSlib tools | ~250MB | latest | amd64, arm64 | 1.22 / 1.22.1 | 1.0.0 |
| ldsc | LD Score Regression for heritability & correlation | ~1.4GB | latest | amd64, arm64 | 2.0.0 | 1.0.0 |

## Versioning Strategy

We use **independent versioning** for our Docker images:
- **Software Version**: The version of the actual software installed in the container
- **Image Version**: Our internal version for the Docker image (starts at 1.0.0)

This allows us to:
- Update the Docker image (bug fixes, optimizations, security updates) while keeping the same software version
- Maintain reproducibility by using specific image versions
- Track changes independently from software releases

Example: `biopsyk/gcta:1.0.0` contains GCTA software version 1.94.1, but if we need to fix the Dockerfile or update dependencies, the next image would be `biopsyk/gcta:1.1.0` with the same GCTA 1.94.1 software.

## Repository Structure

```
.
├── README.md
├── images/
│   ├── plink/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── flashpca2/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── gcta/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── gctb/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── readstat/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── regenie/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   ├── bcftools/
│   │   ├── Dockerfile
│   │   ├── README.md
│   │   └── VERSION
│   └── ldsc/
│       ├── Dockerfile
│       ├── README.md
│       └── VERSION
└── scripts/
    └── build-and-push.sh
```

## Building and Pushing Images

The repository includes a build script that handles multi-architecture builds using Docker BuildX. The script is located in `scripts/build-and-push.sh`.

To build and push an image:

```bash
# Build and push a specific image
./scripts/build-and-push.sh plink

# Build and push all images
./scripts/build-and-push.sh all

# Examples for specific images
./scripts/build-and-push.sh gcta
./scripts/build-and-push.sh gctb
./scripts/build-and-push.sh readstat
./scripts/build-and-push.sh bcftools
./scripts/build-and-push.sh ldsc
```

### Prerequisites

- Docker with BuildX support
- Docker Hub account with access to the institute's organization
- Docker logged in to the registry (`docker login`)

### Building Locally

To build an image locally for your architecture:

```bash
cd images/plink
docker build -t biopsyk/plink:latest .

# Examples for other images
cd images/gcta
docker build -t biopsyk/gcta:latest .

cd images/gctb
docker build -t biopsyk/gctb:latest .

cd images/readstat
docker build -t biopsyk/readstat:latest .

cd images/bcftools
docker build -t biopsyk/bcftools:latest .

cd images/ldsc
docker build -t biopsyk/ldsc:latest .
```

## Contributing

To add a new image:

1. Create a new directory under `images/`
2. Add a Dockerfile and README.md
3. Update the main README.md with image information
4. Test the build locally
5. Submit a pull request

## Using with Singularity

To use these Docker images with Singularity (e.g., on an HPC cluster), you can pull them directly from Docker Hub:

```bash
# Pull using specific version (recommended for reproducibility)
singularity pull docker://biopsyk/plink:1.0.0

# Or pull the latest version (always gets the newest updates)
singularity pull docker://biopsyk/plink:latest

# This will create a Singularity Image File (SIF) named 'plink_1.0.0.sif' or 'plink_latest.sif'
# You can then run PLINK commands using:
singularity exec plink_1.0.0.sif plink1.9 --help
singularity exec plink_1.0.0.sif plink2 --help

# Examples for other tools:
singularity pull docker://biopsyk/gcta:1.94.1
singularity exec gcta_1.94.1.sif gcta64 --help

singularity pull docker://biopsyk/gctb:2.05
singularity exec gctb_2.05.sif gctb --help

singularity pull docker://biopsyk/readstat:1.0.0
singularity exec readstat_1.0.0.sif readstat --help

singularity pull docker://biopsyk/bcftools:1.0.0
singularity exec bcftools_1.0.0.sif bcftools --help
singularity exec bcftools_1.0.0.sif tabix --help

singularity pull docker://biopsyk/ldsc:1.0.0
singularity exec ldsc_1.0.0.sif ldsc -h
singularity exec ldsc_1.0.0.sif munge_sumstats -h

# To bind your data directory (replace /path/to/data with your actual data path):
singularity exec -B /path/to/data:/data plink_1.0.0.sif plink2 --help
```

Note: When using Singularity, the container's `/data` directory is the recommended location for your working files. Use the `-B` flag to bind your host directory to this location.

## License

Please add appropriate license information here.
