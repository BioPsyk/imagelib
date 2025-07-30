# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies (based on working jbn/readstat approach + Ubuntu-specific needs)
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    make \
    cmake \
    git \
    wget \
    zlib1g-dev \
    gettext \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /build

# Build libxlsxwriter first (needed for Excel support)
RUN git clone https://github.com/jmcnamara/libxlsxwriter.git && \
    cd libxlsxwriter && \
    make && \
    make install && \
    cd .. && \
    rm -rf libxlsxwriter

# Clone and build ReadStat (using the working approach from jbn/readstat)
RUN git clone https://github.com/WizardMac/ReadStat.git && \
    cd ReadStat && \
    ./autogen.sh && \
    ./configure && \
    make -j$(nproc) && \
    make install

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    zlib1g \
    && rm -rf /var/lib/apt/lists/*

# Copy the built binaries from builder stage
COPY --from=builder /usr/local/bin/readstat /usr/local/bin/readstat
COPY --from=builder /usr/local/lib/libxlsxwriter* /usr/local/lib/
COPY --from=builder /usr/local/lib/libreadstat* /usr/local/lib/

# Update library cache
RUN ldconfig

# Copy VERSION file and create version information
COPY VERSION /etc/readstat-version-tag
RUN echo "ReadStat CLI - Version $(cat /etc/readstat-version-tag)" > /etc/readstat-version

# Create a simple wrapper script for convenience
RUN echo '#!/bin/bash\n\
# ReadStat Docker wrapper\n\
# Usage: readstat-docker input.dta output.csv\n\
\n\
if [ $# -lt 2 ]; then\n\
    echo "Usage: $0 <input_file> <output_file>"\n\
    echo "Supported input formats: .dta, .por, .sav, .sas7bdat"\n\
    echo "Supported output formats: .dta, .por, .sav, .sas7bdat, .xlsx, .csv"\n\
    exit 1\n\
fi\n\
\n\
readstat "$@"' > /usr/local/bin/readstat-docker && \
    chmod +x /usr/local/bin/readstat-docker

# Set working directory for data
WORKDIR /data

# Default command shows help
CMD ["readstat", "--help"] 