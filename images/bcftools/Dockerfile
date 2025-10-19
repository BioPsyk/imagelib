# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    bzip2 \
    autoconf \
    automake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libgsl-dev \
    libperl-dev \
    && rm -rf /var/lib/apt/lists/*

# Download and build HTSlib
WORKDIR /tmp
RUN curl -L https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 | tar xj && \
    cd htslib-1.22.1 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install

# Download and build bcftools
RUN curl -L https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 | tar xj && \
    cd bcftools-1.22 && \
    ./configure --prefix=/usr/local --enable-libgsl --enable-perl-filters && \
    make -j$(nproc) && \
    make install

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    libcurl4 \
    libssl3 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libncurses5 \
    libgsl27 \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Copy binaries and libraries from builder
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/
COPY --from=builder /usr/local/bin/tabix /usr/local/bin/
COPY --from=builder /usr/local/bin/bgzip /usr/local/bin/
COPY --from=builder /usr/local/bin/htsfile /usr/local/bin/
COPY --from=builder /usr/local/lib/libhts* /usr/local/lib/
COPY --from=builder /usr/local/include/htslib /usr/local/include/htslib
COPY --from=builder /usr/local/share/man/man1/bcftools* /usr/local/share/man/man1/
COPY --from=builder /usr/local/share/man/man1/tabix* /usr/local/share/man/man1/
COPY --from=builder /usr/local/share/man/man1/bgzip* /usr/local/share/man/man1/

# Update library cache
RUN ldconfig

# Copy VERSION file and create version information
COPY VERSION /etc/bcftools-version-tag
RUN echo "BCFtools $(bcftools --version | head -n1 | cut -d' ' -f2) with HTSlib" > /etc/bcftools-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r bcftools && useradd -r -g bcftools bcftools
USER bcftools

# Default command shows help
CMD ["bcftools", "--help"]
