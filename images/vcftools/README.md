# VCFtools Docker Image

This Docker image provides VCFtools (C++ binary and Perl scripts) for working with VCF files.

## Description

VCFtools is a package for processing VCF files such as those produced by the 1000 Genomes Project. It includes the `vcftools` binary plus Perl helpers (`vcf-sort`, `vcf-merge`, `vcf-query`, and others) and the `Vcf.pm` module.

Tabix and bgzip are included because many VCFtools workflows expect bgzip-compressed, tabix-indexed VCF files.

## Version Information

- Internal image version: 1.0.0
- VCFtools version: 0.1.17

## Included Tools

- **vcftools**: C++ binary for filtering, summaries, and related analyses
- **Perl scripts**: `vcf-annotate`, `vcf-compare`, `vcf-concat`, `vcf-consensus`, `vcf-isec`, `vcf-merge`, `vcf-query`, `vcf-sort`, `vcf-stats`, `vcf-subset`, `vcf-to-tab`, `vcf-validator`, and related helpers
- **tabix** / **bgzip**: index and block-compress VCF files

`PERL5LIB` is set so `Vcf.pm` is found automatically.

## Usage

### Basic Usage

```bash
# Summarize a VCF
docker run -v $(pwd):/data biopsyk/vcftools vcftools --vcf input.vcf --out summary

# Filter by minor allele frequency
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --maf 0.01 --recode --out filtered

# Sort a VCF with the Perl helper
docker run -v $(pwd):/data biopsyk/vcftools vcf-sort input.vcf > sorted.vcf

# Compress and index
docker run -v $(pwd):/data biopsyk/vcftools bgzip input.vcf
docker run -v $(pwd):/data biopsyk/vcftools tabix -p vcf input.vcf.gz
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/vcftools:1.0.0

# Run vcftools
singularity exec vcftools_1.0.0.sif vcftools --gzvcf input.vcf.gz --freq --out freqs

# Run a Perl helper
singularity exec vcftools_1.0.0.sif vcf-sort input.vcf > sorted.vcf
```

## Common Analysis Workflows

### 1. Site and genotype summaries

```bash
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --freq --out freqs
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --depth --out depth
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --missing-indv --out missing
```

### 2. Filtering and recoding

```bash
docker run -v $(pwd):/data biopsyk/vcftools vcftools \
  --gzvcf input.vcf.gz \
  --maf 0.01 \
  --max-missing 0.9 \
  --recode --recode-INFO-all \
  --out qc
```

### 3. Relatedness and LD

```bash
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --relatedness2 --out relatedness
docker run -v $(pwd):/data biopsyk/vcftools vcftools --gzvcf input.vcf.gz --hap-r2 --ld-window-bp 100000 --out ld
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- For combined BCFtools + Perl + VCFtools workflows, use the `custom` image
- BCFtools is a separate standalone image (`biopsyk/bcftools`)

## Building the Image

From this directory:

```bash
docker build -t biopsyk/vcftools:latest .
```

## References

- [VCFtools website](https://vcftools.github.io/)
- [VCFtools GitHub repository](https://github.com/vcftools/vcftools)
- [Usage examples](https://vcftools.github.io/examples.html)
- Danecek P, et al. (2011) The variant call format and VCFtools. Bioinformatics, 27(15):2156-2158.
