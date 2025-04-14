# Stage 1: Build environment
FROM ubuntu:22.04 as builder

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    zlib1g-dev \
    libopenblas-dev \
    libatlas-base-dev \
    cmake \
    libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

# Create and set working directory
WORKDIR /plink

# Download and extract source code
RUN wget https://github.com/chrchang/plink-ng/archive/refs/heads/master.zip && \
    unzip master.zip && \
    rm master.zip && \
    ls -la plink-ng-master/2.0

# Build PLINK 1.9
WORKDIR /plink/plink-ng-master/1.9
RUN sed -i 's/^#define NOLAPACK/#define USE_OPENBLAS\n\/\/#define NOLAPACK/' plink_common.h && \
    sed -i 's/^#define NO_SSE2/\/\/#define NO_SSE2/' plink_common.h && \
    make FORCE_DYNAMIC=1 ZLIB=-lz BLASFLAGS="-lopenblas" CFLAGS='-Wall -O2 -I. -I../2.0/simde -fPIC' CXXFLAGS='-Wall -O2 -I. -I../2.0/simde -fPIC' plink && \
    cp plink /usr/local/bin/plink1.9 && \
    plink1.9 --help | head -n 1

# Build PLINK 2.0
WORKDIR /plink/plink-ng-master/2.0
RUN ls -la && \
    make -j$(nproc) ZLIB=-lz && \
    cp bin/plink2 /usr/local/bin/plink2

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    zlib1g \
    libopenblas0 \
    libatlas3-base \
    && rm -rf /var/lib/apt/lists/*

# Copy the built binary from builder stage
COPY --from=builder /usr/local/bin/plink1.9 /usr/local/bin/plink1.9
COPY --from=builder /usr/local/bin/plink2 /usr/local/bin/plink2

# Create version file
RUN echo "PLINK 1.9 and PLINK 2.0 (built from source)" > /etc/plink-version

# Set working directory for data
WORKDIR /data

# Verify installation
RUN plink1.9 --help | head -n 1

# Verify installations and provide version info
RUN plink2 --version >> /etc/plink-version

# Default command
CMD ["plink1.9", "--help"] 