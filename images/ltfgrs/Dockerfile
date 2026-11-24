# LTFGRS: phenotype-based family genetic risk scores (LT-FH++, PA-FGRS, Kendler FGRS)
# Upstream: https://github.com/EmilMiP/LTFGRS (latest main)
# Docs: https://emilmip.github.io/LTFGRS/

FROM rocker/r-ver:4.4.2

ARG DEBIAN_FRONTEND=noninteractive
ARG LTFGRS_REF=main

LABEL org.opencontainers.image.title="LTFGRS" \
      org.opencontainers.image.description="Phenotype-based family genetic risk scores (LT-FH++, PA-FGRS, Kendler FGRS)" \
      org.opencontainers.image.source="https://github.com/EmilMiP/LTFGRS"

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        build-essential \
        cmake \
        gfortran \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libglpk-dev \
        libgmp-dev \
        libblas-dev \
        liblapack-dev \
        libicu-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'

WORKDIR /tmp
RUN git clone --depth 1 --branch "${LTFGRS_REF}" https://github.com/EmilMiP/LTFGRS.git LTFGRS \
    && mkdir -p /opt/ltfgrs \
    && git -C /tmp/LTFGRS rev-parse HEAD > /opt/ltfgrs/IMAGE_REVISION \
    && Rscript -e 'options(Ncpus = max(1L, parallel::detectCores())); remotes::install_local("/tmp/LTFGRS", dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never")' \
    && Rscript -e 'install.packages("kinship2", repos = "https://cloud.r-project.org")' \
    && rm -rf /tmp/LTFGRS /tmp/Rtmp*

# Keep upstream vignettes as plain .Rmd (no pandoc/rmarkdown/HTML).
RUN git clone --depth 1 --branch "${LTFGRS_REF}" https://github.com/EmilMiP/LTFGRS.git /tmp/LTFGRS \
    && cp -a /tmp/LTFGRS/vignettes /opt/ltfgrs/vignettes \
    && rm -rf /tmp/LTFGRS

COPY VERSION /etc/ltfgrs-version-tag
RUN { \
      echo "LTFGRS ${LTFGRS_REF} ($(tr -d '[:space:]' < /opt/ltfgrs/IMAGE_REVISION))"; \
      echo "Package $(Rscript -e "cat(as.character(packageVersion('LTFGRS')))")"; \
      echo "Image $(cat /etc/ltfgrs-version-tag)"; \
      Rscript -e 'cat(R.version.string, "\n")'; \
    } > /etc/ltfgrs-version \
    && printf '%s\n' '#!/bin/sh' 'exec cat /etc/ltfgrs-version' > /usr/local/bin/ltfgrs-version \
    && chmod 0755 /usr/local/bin/ltfgrs-version

ENV HOME=/tmp/ltfgrs-home

WORKDIR /data

RUN groupadd --system ltfgrs \
    && useradd --system --gid ltfgrs --home-dir /tmp/ltfgrs-home ltfgrs \
    && mkdir -p /tmp/ltfgrs-home /data \
    && chown -R ltfgrs:ltfgrs /tmp/ltfgrs-home /data /opt/ltfgrs
USER ltfgrs

CMD ["ltfgrs-version"]
