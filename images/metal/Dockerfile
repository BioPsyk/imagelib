# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    ca-certificates \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone METAL repository
WORKDIR /tmp
RUN git clone https://github.com/statgen/METAL.git && \
    cd METAL && \
    git checkout $(git describe --tags --abbrev=0 2>/dev/null || echo "master")

# Build METAL
WORKDIR /tmp/METAL
RUN mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make && \
    make test

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy METAL binary from builder
COPY --from=builder /tmp/METAL/build/metal/metal /usr/local/bin/metal

# Copy VERSION file and create version information
COPY VERSION /etc/metal-version-tag
RUN metal --version 2>&1 | head -1 > /etc/metal-version 2>/dev/null || \
    echo "METAL Meta-Analysis Tool - Version $(cat /etc/metal-version-tag)" > /etc/metal-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r metal && useradd -r -g metal metal
USER metal

# Default command shows version and help
CMD ["metal", "--help"]
