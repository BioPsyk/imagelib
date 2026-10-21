# GetPheno: extract phenotypes from biobank-scale files
# Upstream: https://github.com/MischaLundberg/GetPheno (V1.30)

FROM rocker/r-ver:4.3.3

ARG DEBIAN_FRONTEND=noninteractive
ARG GETPHENO_REF=V1.30
ARG GETPHENO_SHA=f49d73ae1c3a5868d7d2470d27ad61315993951a

LABEL org.opencontainers.image.source="https://github.com/MischaLundberg/GetPheno" \
      org.opencontainers.image.revision="${GETPHENO_SHA}"

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
    && test "$(git rev-parse HEAD)" = "${GETPHENO_SHA}" \
    && python3 -m pip install --no-cache-dir -r docker/requirements.txt \
    && Rscript -e 'install.packages("shiny", repos = getOption("repos"))' \
    && mkdir -p /opt \
    && mv /tmp/getpheno /opt/getpheno \
    && chmod 0755 /opt/getpheno/get_pheno.py /opt/getpheno/docker/entrypoint.sh \
    && python3 -m py_compile /opt/getpheno/get_pheno.py \
    && printf '%s\n' "${GETPHENO_SHA}" > /opt/getpheno/IMAGE_REVISION \
    && printf '#!/bin/bash\nexec python3 /opt/getpheno/get_pheno.py "$@"\n' > /usr/local/bin/get_pheno.py \
    && printf '#!/bin/bash\nexec python3 /opt/getpheno/get_pheno.py "$@"\n' > /usr/local/bin/getpheno \
    && printf '#!/bin/bash\nexec /opt/getpheno/docker/entrypoint.sh "$@"\n' > /usr/local/bin/getpheno-shiny \
    && chmod 0755 /usr/local/bin/get_pheno.py /usr/local/bin/getpheno /usr/local/bin/getpheno-shiny

ENV GETPHENO_SHINY_CONFIG=/config/projects.conf \
    GETPHENO_HOST=0.0.0.0 \
    GETPHENO_PORT=3838 \
    HOME=/tmp/getpheno-home \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY VERSION /etc/getpheno-version-tag
RUN { \
      echo "GetPheno ${GETPHENO_REF} (${GETPHENO_SHA})"; \
      echo "Image $(cat /etc/getpheno-version-tag)"; \
      python3 --version; \
      Rscript -e 'cat(R.version.string, "\n")'; \
    } > /etc/getpheno-version

WORKDIR /data

RUN groupadd -r getpheno && useradd -r -g getpheno getpheno \
    && mkdir -p /tmp/getpheno-home \
    && chown -R getpheno:getpheno /tmp/getpheno-home /opt/getpheno /data
USER getpheno

EXPOSE 3838

CMD ["getpheno", "--help"]
