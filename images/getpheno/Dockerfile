# GetPheno: extract phenotypes from biobank-scale files
# Upstream: https://github.com/MischaLundberg/GetPheno (main)

FROM rocker/r-ver:4.3.3

ARG DEBIAN_FRONTEND=noninteractive
ARG GETPHENO_REF=main

LABEL org.opencontainers.image.title="GetPheno" \
      org.opencontainers.image.description="Extract phenotypes from biobank-scale files" \
      org.opencontainers.image.source="https://github.com/MischaLundberg/GetPheno"

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        python3 \
        python3-pip \
        python3-dev \
        build-essential \
        pkg-config \
        libhdf5-dev \
        liblzo2-dev \
        libbz2-dev \
        slurm-client \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN git clone --depth 1 --branch "${GETPHENO_REF}" https://github.com/MischaLundberg/GetPheno.git getpheno \
    && cd getpheno \
    && python3 -m pip install --no-cache-dir -r docker/requirements.txt \
    && Rscript -e 'install.packages("shiny", repos = getOption("repos"))' \
    && mkdir -p /opt \
    && git rev-parse HEAD > /tmp/getpheno-revision \
    && mv /tmp/getpheno /opt/getpheno \
    && mv /tmp/getpheno-revision /opt/getpheno/IMAGE_REVISION \
    && chmod 0755 /opt/getpheno/get_pheno.py /opt/getpheno/docker/entrypoint.sh \
    && python3 -m py_compile /opt/getpheno/get_pheno.py \
    && printf '#!/bin/sh\nexec python3 /opt/getpheno/get_pheno.py "$@"\n' > /usr/local/bin/get_pheno.py \
    && printf '#!/bin/sh\nexec python3 /opt/getpheno/get_pheno.py "$@"\n' > /usr/local/bin/getpheno \
    && printf '#!/bin/sh\nexec /opt/getpheno/docker/entrypoint.sh shiny "$@"\n' > /usr/local/bin/getpheno-shiny \
    && chmod 0755 /usr/local/bin/get_pheno.py /usr/local/bin/getpheno /usr/local/bin/getpheno-shiny

ENV GETPHENO_SHINY_CONFIG=/config/projects.conf \
    GETPHENO_HOST=0.0.0.0 \
    GETPHENO_PORT=3838 \
    HOME=/tmp/getpheno-home \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY VERSION /etc/getpheno-version-tag
RUN { \
      echo "GetPheno ${GETPHENO_REF} ($(tr -d '[:space:]' < /opt/getpheno/IMAGE_REVISION))"; \
      echo "Image $(cat /etc/getpheno-version-tag)"; \
      python3 --version; \
      Rscript -e 'cat(R.version.string, "\n")'; \
    } > /etc/getpheno-version

WORKDIR /data

RUN groupadd --system getpheno \
    && useradd --system --gid getpheno --home-dir /tmp/getpheno-home getpheno \
    && mkdir -p /tmp/getpheno-home /data \
    && chown -R getpheno:getpheno /tmp/getpheno-home /data /opt/getpheno
USER getpheno

EXPOSE 3838

# CLI-first: upstream defaults to Shiny. Their entrypoint still accepts getpheno/--help.
ENTRYPOINT ["/opt/getpheno/docker/entrypoint.sh"]
CMD ["getpheno", "--help"]
