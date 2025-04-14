# Institute Docker Images

This repository contains Dockerfiles for commonly used software at our institute. Each image is optimized for size and built for multiple architectures (amd64 and arm64).

## Available Images

| Image Name | Description | Size | Latest Tag | Architectures |
|------------|-------------|------|------------|---------------|
| plink | PLINK1.9 and PLINK2 for genetic analysis | ~100MB | latest | amd64, arm64 |

## Repository Structure

```
.
├── README.md
├── images/
│   └── plink/
│       ├── Dockerfile
│       └── README.md
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
# Pull the latest PLINK image
singularity pull docker://biopsyk/plink:latest

# This will create a Singularity Image File (SIF) named 'plink_latest.sif'
# You can then run PLINK commands using:
singularity exec plink_latest.sif plink1.9 --help
singularity exec plink_latest.sif plink2 --help

# To bind your data directory (replace /path/to/data with your actual data path):
singularity exec -B /path/to/data:/data plink_latest.sif plink2 --help
```

Note: When using Singularity, the container's `/data` directory is the recommended location for your working files. Use the `-B` flag to bind your host directory to this location.

## License

Please add appropriate license information here.
