# pedsum: pedigree summary CLI
# Upstream: https://github.com/rwaples/pedsum (latest main)
# Requires Python 3.14; 3.15 is not supported yet (polars 1.x string methods).

FROM python:3.14-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive
ARG PEDSUM_REF=main

LABEL org.opencontainers.image.title="pedsum" \
      org.opencontainers.image.description="Pedigree summary CLI: structure, relatedness, inbreeding, and Ne estimators" \
      org.opencontainers.image.source="https://github.com/rwaples/pedsum"

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        pigz \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp
RUN git clone --depth 1 --branch "${PEDSUM_REF}" https://github.com/rwaples/pedsum.git pedsum \
    && python -m pip install --no-cache-dir --upgrade pip \
    && python -m pip install --no-cache-dir ./pedsum \
    && mkdir -p /opt/pedsum \
    && cp /tmp/pedsum/example_pedigree.tsv /opt/pedsum/example_pedigree.tsv \
    && git -C /tmp/pedsum rev-parse HEAD > /opt/pedsum/IMAGE_REVISION \
    && rm -rf /tmp/pedsum

COPY VERSION /etc/pedsum-version-tag
RUN { \
      echo "pedsum ${PEDSUM_REF} ($(tr -d '[:space:]' < /opt/pedsum/IMAGE_REVISION))"; \
      echo "Image $(cat /etc/pedsum-version-tag)"; \
      python --version; \
      pedsum --version; \
    } > /etc/pedsum-version

WORKDIR /data

RUN groupadd --system pedsum \
    && useradd --system --gid pedsum --home-dir /tmp/pedsum-home pedsum \
    && mkdir -p /tmp/pedsum-home /data \
    && chown -R pedsum:pedsum /tmp/pedsum-home /data /opt/pedsum
USER pedsum

CMD ["pedsum", "--help"]
