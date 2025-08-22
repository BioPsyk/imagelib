# Use Ubuntu as base image - x86_64 only
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    unzip \
    libgomp1 \
    libstdc++6 \
    libgcc-s1 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libboost-iostreams1.74.0 \
    || true && \
    rm -rf /var/lib/apt/lists/*

# Download REGENIE binary directly from GitHub releases
# Using version 3.6 which is more stable and widely tested
RUN curl -kL https://github.com/rgcgithub/regenie/releases/download/v3.6/regenie_v3.6.gz_x86_64_Linux_mkl.zip -o /tmp/regenie.zip && \
    cd /tmp && \
    unzip -q regenie.zip && \
    mv regenie_v3.6.gz_x86_64_Linux_mkl /usr/local/bin/regenie && \
    chmod +x /usr/local/bin/regenie && \
    rm /tmp/regenie.zip

# Copy VERSION file and create version information
COPY VERSION /etc/regenie-version-tag
RUN echo "REGENIE (REGression with ENvironment Interaction and polygenic Effects) - Version 3.6" > /etc/regenie-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r regenie && useradd -r -g regenie regenie
USER regenie

# Default command shows help
CMD ["regenie", "--help"]