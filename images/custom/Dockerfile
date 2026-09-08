# Combined image: BCFtools + HTSlib tools, Perl, and VCFtools
# Software versions match the standalone images.

# Stage 1: Build HTSlib and BCFtools
FROM ubuntu:22.04 AS bcftools-builder

ENV DEBIAN_FRONTEND=noninteractive

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

WORKDIR /tmp
RUN curl -L https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 | tar xj && \
    cd htslib-1.22.1 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install

RUN curl -L https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 | tar xj && \
    cd bcftools-1.22 && \
    ./configure --prefix=/usr/local --enable-libgsl --enable-perl-filters && \
    make -j$(nproc) && \
    make install

# Stage 2: Build VCFtools
FROM ubuntu:22.04 AS vcftools-builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    pkg-config \
    zlib1g-dev \
    perl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN curl -L https://github.com/vcftools/vcftools/releases/download/v0.1.17/vcftools-0.1.17.tar.gz | tar xz && \
    cd vcftools-0.1.17 && \
    ./configure --prefix=/opt/vcftools --with-pmdir=share/perl5 && \
    make -j$(nproc) && \
    make install

# Stage 3: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PERL5LIB=/usr/local/share/perl5

# Runtime libraries for BCFtools/HTSlib plus a full Perl for VCFtools scripts
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

# BCFtools and HTSlib utilities
COPY --from=bcftools-builder /usr/local/bin/bcftools /usr/local/bin/
COPY --from=bcftools-builder /usr/local/bin/tabix /usr/local/bin/
COPY --from=bcftools-builder /usr/local/bin/bgzip /usr/local/bin/
COPY --from=bcftools-builder /usr/local/bin/htsfile /usr/local/bin/
COPY --from=bcftools-builder /usr/local/lib/libhts* /usr/local/lib/
COPY --from=bcftools-builder /usr/local/include/htslib /usr/local/include/htslib

# VCFtools C++ binary, Perl scripts, and Vcf.pm module
COPY --from=vcftools-builder /opt/vcftools/bin/ /usr/local/bin/
COPY --from=vcftools-builder /opt/vcftools/share/perl5 /usr/local/share/perl5

RUN ldconfig

COPY VERSION /etc/custom-version-tag
RUN { \
      echo "custom image $(cat /etc/custom-version-tag)"; \
      echo "BCFtools $(bcftools --version | head -n1 | cut -d' ' -f2) with HTSlib"; \
      echo "Perl $(perl -e 'print $^V')"; \
      vcftools --version; \
    } > /etc/custom-version

WORKDIR /data

RUN groupadd -r custom && useradd -r -g custom custom
USER custom

CMD ["bcftools", "--help"]
