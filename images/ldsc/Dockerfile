# Stage 1: Build environment
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    git \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 - && \
    ln -s /root/.local/bin/poetry /usr/local/bin/poetry

# Clone LDSC Python 3 version
WORKDIR /tmp
RUN git clone https://github.com/abdenlab/ldsc-python3.git && \
    cd ldsc-python3 && \
    git checkout $(git describe --tags --abbrev=0 2>/dev/null || echo "main")

# Install LDSC dependencies
WORKDIR /tmp/ldsc-python3
ENV POETRY_VIRTUALENVS_CREATE=false
RUN poetry install --only main --no-interaction --no-root && \
    poetry install --only main --no-interaction

# Stage 2: Runtime environment
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Copy Python packages from builder
COPY --from=builder /usr/local/lib/python3.10/dist-packages /usr/local/lib/python3.10/dist-packages

# Copy LDSC application
COPY --from=builder /tmp/ldsc-python3 /usr/local/lib/ldsc
WORKDIR /usr/local/lib/ldsc

# Create wrapper scripts in PATH
RUN echo '#!/bin/bash\npython3 /usr/local/lib/ldsc/ldsc.py "$@"' > /usr/local/bin/ldsc && \
    echo '#!/bin/bash\npython3 /usr/local/lib/ldsc/munge_sumstats.py "$@"' > /usr/local/bin/munge_sumstats && \
    echo '#!/bin/bash\npython3 /usr/local/lib/ldsc/make_annot.py "$@"' > /usr/local/bin/make_annot && \
    chmod +x /usr/local/bin/ldsc /usr/local/bin/munge_sumstats /usr/local/bin/make_annot

# Copy VERSION file and create version information
COPY VERSION /etc/ldsc-version-tag
RUN python3 -c "import ldscore; print('LDSC (LD Score Regression) - Version:', open('/etc/ldsc-version-tag').read().strip())" > /etc/ldsc-version 2>/dev/null || \
    echo "LDSC (LD Score Regression) - Version $(cat /etc/ldsc-version-tag)" > /etc/ldsc-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r ldsc && useradd -r -g ldsc ldsc && \
    chown -R ldsc:ldsc /usr/local/lib/ldsc
USER ldsc

# Default command shows help
CMD ["ldsc", "-h"]
