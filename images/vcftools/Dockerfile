# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    pkg-config \
    zlib1g-dev \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Download and build VCFtools (C++ binary + Perl scripts/modules)
WORKDIR /tmp
RUN curl -L https://github.com/vcftools/vcftools/releases/download/v0.1.17/vcftools-0.1.17.tar.gz | tar xz && \
    cd vcftools-0.1.17 && \
    ./configure --prefix=/opt/vcftools --with-pmdir=share/perl5 && \
    make -j$(nproc) && \
    make install

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PERL5LIB=/usr/local/share/perl5

# Perl (for vcf-* scripts), zlib, and tabix/bgzip used by many VCFtools workflows
RUN apt-get update && apt-get install -y \
    perl \
    zlib1g \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy VCFtools binary, Perl scripts, modules, and man pages
COPY --from=builder /opt/vcftools/bin/ /usr/local/bin/
COPY --from=builder /opt/vcftools/share/perl5 /usr/local/share/perl5
COPY --from=builder /opt/vcftools/share/man /usr/local/share/man

# Copy VERSION file and create version information
COPY VERSION /etc/vcftools-version-tag
RUN vcftools --version > /etc/vcftools-version && \
    echo "Image $(cat /etc/vcftools-version-tag)" >> /etc/vcftools-version

WORKDIR /data

# Create a non-root user for security
RUN groupadd -r vcftools && useradd -r -g vcftools vcftools
USER vcftools

CMD ["vcftools", "--help"]
