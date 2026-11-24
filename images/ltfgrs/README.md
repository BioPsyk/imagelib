# LTFGRS Docker Image

This Docker image provides [LTFGRS](https://github.com/EmilMiP/LTFGRS), an R package for phenotype-based family genetic risk scores: LT-FH++, PA-FGRS, and the simplified Kendler FGRS. Documentation: https://emilmip.github.io/LTFGRS/

## Description

LTFGRS estimates genetic liability from family history using personalised age-, sex-, and birth-year-specific thresholds. Core functions include `prepare_thresholds`, `estimate_liability`, `prepare_graph`, `familywise_censoring`, and `kendler_simplified`. Parallel estimation uses the `future` package.

The image is an R library, not a CLI. Run `R` or `Rscript` and `library(LTFGRS)`.

## Version Information

- Internal image version: 1.1.0
- LTFGRS version: latest `main` at image build time (GitHub 1.2.0; CRAN is 1.0.1). Exact commit is written to `/opt/ltfgrs/IMAGE_REVISION` and `/etc/ltfgrs-version`
- R 4.4.2 (rocker/r-ver)
- Extra Suggests package: `kinship2` (used in the upstream tutorials)

Installed versions are also recorded in `/etc/ltfgrs-version`.

## Included Tools

- **R** / **Rscript**: run analyses that `library(LTFGRS)`
- **ltfgrs-version**: print the image and package versions
- LTFGRS Imports (including Rcpp, igraph, xgboost, future, ggplot2, dplyr)
- Upstream vignette sources as `.Rmd` under `/opt/ltfgrs/vignettes` (plain text; not knitted HTML)

## Usage

### Basic Usage

```bash
# Show installed versions
docker run --rm biopsyk/ltfgrs

# Interactive R
docker run --rm -it -v $(pwd):/data biopsyk/ltfgrs R

# Run a script (mount input/output into /data)
docker run --rm -v $(pwd):/data biopsyk/ltfgrs Rscript analysis.R
```

Example session:

```r
library(LTFGRS)
packageVersion("LTFGRS")
# prepare_thresholds(), estimate_liability(), prepare_graph(), kendler_simplified()
```

For parallel runs inside the container:

```r
library(future)
plan(multisession)
```

Rendered HTML tutorials: https://emilmip.github.io/LTFGRS/articles/

### Vignette sources (plain `.Rmd`)

The container does not build HTML vignettes (no pandoc/rmarkdown). The upstream `.Rmd` files are in `/opt/ltfgrs/vignettes` and can be copied out and opened in an editor, or sourced as R:

```bash
# List them
docker run --rm biopsyk/ltfgrs ls /opt/ltfgrs/vignettes

# Copy into the current directory
docker run --rm -v $(pwd):/data biopsyk/ltfgrs \
  cp -a /opt/ltfgrs/vignettes /data/ltfgrs-vignettes
```

### Using with Singularity

```bash
singularity pull docker://biopsyk/ltfgrs:1.1.0
singularity exec ltfgrs_1.1.0.sif ltfgrs-version
singularity exec ltfgrs_1.1.0.sif ls /opt/ltfgrs/vignettes
singularity exec -B /path/to/data:/data ltfgrs_1.1.0.sif \
  cp -a /opt/ltfgrs/vignettes /data/ltfgrs-vignettes
singularity exec -B /path/to/data:/data ltfgrs_1.1.0.sif Rscript /data/analysis.R
singularity run -B /path/to/data:/data ltfgrs_1.1.0.sif
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- Rebuild to pick up a newer GitHub `main`; the cloned commit is in `/opt/ltfgrs/IMAGE_REVISION`

## Building the Image

From this directory:

```bash
docker build -t biopsyk/ltfgrs:latest .
```

To pin a different upstream ref:

```bash
docker build --build-arg LTFGRS_REF=<tag-or-branch> -t biopsyk/ltfgrs:latest .
```
