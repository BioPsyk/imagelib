# pedsum Docker Image

This Docker image provides [pedsum](https://github.com/rwaples/pedsum), a Python CLI for summarizing and validating pedigrees.

## Description

pedsum reports pedigree size and structure, sibship-size distribution, relationship-pair counts (up to degree 5), mating-pair structure, founder contributions, per-individual genealogy, and estimates effective population size eight ways. It also validates pedigree integrity, estimates mate correlation, and can emit EPIMIGHT input.

Relationship-pair detection uses [`pedigree-graph`](https://github.com/rwaples/pedigree-graph). Output terminology follows the upstream [GLOSSARY.md](https://github.com/rwaples/pedsum/blob/main/GLOSSARY.md).

## Version Information

- Internal image version: 1.1.0
- pedsum version: latest `main` at image build time (exact commit is written to `/opt/pedsum/IMAGE_REVISION` and `/etc/pedsum-version`)
- Python 3.14 (required; 3.15 is not supported yet)

Installed versions are also recorded in `/etc/pedsum-version`. The cloned commit is written to `/opt/pedsum/IMAGE_REVISION`.

## Included Tools

- **pedsum**: pedigree summary CLI (`validate`, `summarize`, `effective-size`, `assortative-mating`, `epimight-input`)
- **pigz**: faster `.tsv.gz` writes (pedsum falls back to gzip if missing)
- Example pedigree: `/opt/pedsum/example_pedigree.tsv`

## Usage

### Basic Usage

```bash
# Show help
docker run -v $(pwd):/data biopsyk/pedsum pedsum --help

# Validate a pedigree
docker run -v $(pwd):/data biopsyk/pedsum pedsum validate --in pedigree.tsv --out validate_out

# Summarize (writes summary.yaml, summary.extra.yaml, annotated.tsv.gz)
docker run -v $(pwd):/data biopsyk/pedsum pedsum summarize --in pedigree.tsv --out summary_out

# Effective population size
docker run -v $(pwd):/data biopsyk/pedsum pedsum effective-size \
  --in pedigree.tsv --out ne_out --birth-year-col birth_year

# Mate correlation
docker run -v $(pwd):/data biopsyk/pedsum pedsum assortative-mating \
  --in pedigree.tsv --out am_out --trait liability

# EPIMIGHT input
docker run -v $(pwd):/data biopsyk/pedsum pedsum epimight-input \
  --in pedigree.tsv --out epimight_out --pairs
```

`--out` is a directory (created if needed). Required pedigree columns are `id`, `sex`, `mother`, and `father` (PLINK numeric sex: `1` = male, `2` = female). Extra columns are preserved in `annotated.tsv.gz`.

### Example pedigree

The image includes the 200-individual, 5-generation example from the pedsum repo:

```bash
docker run --rm biopsyk/pedsum pedsum validate \
  --in /opt/pedsum/example_pedigree.tsv --out /tmp/demo-validate
```

### Using with Singularity

```bash
singularity pull docker://biopsyk/pedsum:1.1.0
singularity exec pedsum_1.1.0.sif pedsum --help
singularity exec -B /path/to/data:/data pedsum_1.1.0.sif \
  pedsum summarize --in /data/pedigree.tsv --out /data/summary_out
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- `summarize` computes inbreeding (F) by default; pass `--no-inbreeding` to skip it
- Run `validate` first if the input may have cycles, missing parents, or sex conflicts
- Full command reference: https://github.com/rwaples/pedsum

## Building the Image

From this directory:

```bash
docker build -t biopsyk/pedsum:latest .
```

To pin a different upstream tag:

```bash
docker build --build-arg PEDSUM_REF=main -t biopsyk/pedsum:latest .
```

## References

- [pedsum GitHub repository](https://github.com/rwaples/pedsum)
- [pedigree-graph](https://github.com/rwaples/pedigree-graph)
- [pedsum glossary](https://github.com/rwaples/pedsum/blob/main/GLOSSARY.md)

## License

pedsum is MIT-licensed. See the LICENSE file in the GitHub repository.
