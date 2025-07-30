# Stage 1: Build environment with conda
# Use x86_64 conda regardless of target platform since GCTB is only available for x86_64
FROM --platform=linux/amd64 continuumio/miniconda3 AS builder

# Install GCTB from bioconda (pre-built binary for x86_64)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y gctb && \
    conda clean --all -y

# Stage 2: Runtime environment (respects target platform)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the GCTB binary from conda environment
COPY --from=builder /opt/conda/bin/gctb /usr/local/bin/gctb

# Copy VERSION file and create version information
COPY VERSION /etc/gctb-version-tag
RUN echo "GCTB (Genome-wide Complex Trait Bayesian analysis) - Version $(cat /etc/gctb-version-tag)" > /etc/gctb-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r gctb && useradd -r -g gctb gctb
USER gctb

# Default command shows help
CMD ["gctb", "--help"] 