# Use Ubuntu 22.04 as base image
FROM ubuntu:22.04 AS builder

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    make \
    git \
    wget \
    libmd0 \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libeigen3-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Download and extract Spectra library
WORKDIR /opt
RUN wget https://github.com/yixuan/spectra/archive/v0.8.1.tar.gz && \
    tar xvf v0.8.1.tar.gz && \
    rm v0.8.1.tar.gz

# Clone and build flashpca
WORKDIR /opt
RUN git clone https://github.com/gabraham/flashpca.git && \
    cd flashpca && \
    make all \
    EIGEN_INC=/usr/include/eigen3 \
    BOOST_INC=/usr/include/boost \
    SPECTRA_INC=/opt/spectra-0.8.1/include

# Create final minimal image
FROM ubuntu:22.04

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libboost-program-options1.74.0 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled binary from builder
COPY --from=builder /opt/flashpca/flashpca /usr/local/bin/

# Copy VERSION file and create version information
COPY VERSION /etc/flashpca-version-tag
RUN echo "flashpca2 (built from source) - Version $(cat /etc/flashpca-version-tag)" > /etc/flashpca-version

# Create a non-root user
RUN useradd -m flashpca && \
    mkdir -p /data && \
    chown flashpca:flashpca /data

# Set working directory and user
WORKDIR /data
USER flashpca

# Set the entrypoint
ENTRYPOINT ["flashpca"] 