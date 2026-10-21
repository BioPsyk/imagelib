# GetPheno Docker Image

This Docker image provides [GetPheno](https://github.com/MischaLundberg/GetPheno), a Python tool for extracting phenotypes from large, biobank-scale files (DST, Danish registers, DBDS, CHB, iPSYCH, BioVU, ATLAS, UKBB, and similar).

## Description

GetPheno builds one or more case/control phenotypes from ICD8/9/10 (including CM) and ATC code requests. It supports multi-phenotype request files, advanced `main=` / `sub=` / `rule_out=` definitions, exclusion files, IID subsets, low-memory batching, and optional OMOP CDM cohort output.

The image also includes the optional R Shiny interface and a Slurm client, matching the upstream Docker image. The default command is the CLI (`getpheno --help`), not the Shiny server.

## Version Information

- Internal image version: 1.0.0
- GetPheno version: 1.30 (`V1.30`, commit `f49d73ae1c3a5868d7d2470d27ad61315993951a`)
- Python 3 and R 4.3.3 (rocker/r-ver)

Installed versions are also recorded in `/etc/getpheno-version`.

## Included Tools

- **getpheno** / **get_pheno.py**: phenotype builder CLI
- **getpheno-shiny**: optional Shiny UI (requires a mounted project registry)
- Slurm client tools for cluster job submission from the Shiny app

## Usage

### Basic Usage

```bash
# Show CLI help
docker run -v $(pwd):/data biopsyk/getpheno getpheno --help

# Build phenotypes (mount input/output into /data)
docker run -v $(pwd):/data biopsyk/getpheno getpheno \
  -g diagnosis.request \
  -o diagnosis.casecontrol

# Optional INI in the working directory overrides defaults; CLI flags win
docker run -v $(pwd):/data biopsyk/getpheno getpheno \
  -g phenotypes.tsv \
  -f diagnoses.tsv \
  -i cohort.tsv \
  -o output.tsv
```

### Shiny app

The upstream entrypoint starts Shiny and expects a readable project registry at `/config/projects.conf`. See the [GetPheno Docker deployment docs](https://github.com/MischaLundberg/GetPheno/blob/main/shiny_app/README.md#docker-deployment).

```bash
docker run --rm -p 3838:3838 \
  -v /path/to/projects.conf:/config/projects.conf:ro \
  biopsyk/getpheno getpheno-shiny
```

The official compose-based cluster deployment remains `ghcr.io/mischalundberg/getpheno`.

### Using with Singularity

```bash
singularity pull docker://biopsyk/getpheno:1.0.0
singularity exec getpheno_1.0.0.sif getpheno --help
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- Default matching is prefix-based; use `--eM` for exact matches except codes with a trailing `*`
- `get_pheno.ini` is read from the current working directory or from `/opt/getpheno`

## Building the Image

From this directory:

```bash
docker build -t biopsyk/getpheno:latest .
```

## References

- [GetPheno GitHub repository](https://github.com/MischaLundberg/GetPheno)
- [GetPheno V1.30](https://github.com/MischaLundberg/GetPheno/releases/tag/V1.30)
- [DOI 10.5281/zenodo.21902133](https://doi.org/10.5281/zenodo.21902133)
- Upstream image: `ghcr.io/mischalundberg/getpheno`
