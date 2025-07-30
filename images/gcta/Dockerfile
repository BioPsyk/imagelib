# Stage 1: Build environment with conda
# Use x86_64 conda regardless of target platform since GCTA is only available for x86_64
FROM --platform=linux/amd64 continuumio/miniconda3 AS builder

# Install GCTA from bioconda (pre-built binary for x86_64)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y gcta && \
    conda clean --all -y

# Stage 2: Runtime environment (respects target platform)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the GCTA binary from conda environment
COPY --from=builder /opt/conda/bin/gcta64 /usr/local/bin/gcta64

# Create symlink for convenience
RUN ln -s /usr/local/bin/gcta64 /usr/local/bin/gcta

# Copy VERSION file and create version information
COPY VERSION /etc/gcta-version-tag
RUN echo "GCTA (Genome-wide Complex Trait Analysis) - Version $(cat /etc/gcta-version-tag)" > /etc/gcta-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r gcta && useradd -r -g gcta gcta
USER gcta

# Default command shows help
CMD ["gcta64", "--help"] 