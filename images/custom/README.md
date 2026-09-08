# Custom Docker Image (BCFtools + Perl + VCFtools)

This image combines the three tools that are also published as standalone imagelib images:

- `biopsyk/bcftools`
- `biopsyk/perl`
- `biopsyk/vcftools`

Use it when a single container must run BCFtools, Perl scripts, and VCFtools in one workflow.

## Description

The custom image is a convenience composition, not a different software stack. Software versions match the standalone images:

- **BCFtools 1.22** with HTSlib 1.22.1 (`tabix`, `bgzip`, `htsfile`)
- **Perl 5.34** (Ubuntu 22.04)
- **VCFtools 0.1.17** (C++ binary, Perl helpers, and `Vcf.pm`)

`PERL5LIB` includes the VCFtools Perl modules so helpers such as `vcf-sort` and `vcf-merge` work without extra setup.

## Version Information

- Internal image version: 1.0.0
- BCFtools version: 1.22
- HTSlib version: 1.22.1
- Perl version: 5.34
- VCFtools version: 0.1.17

Installed versions are also recorded in `/etc/custom-version`.

## Usage

### Basic Usage

```bash
# BCFtools
docker run -v $(pwd):/data biopsyk/custom bcftools view input.vcf.gz

# VCFtools
docker run -v $(pwd):/data biopsyk/custom vcftools --gzvcf input.vcf.gz --freq --out freqs

# Perl
docker run -v $(pwd):/data biopsyk/custom perl myscript.pl

# VCFtools Perl helper
docker run -v $(pwd):/data biopsyk/custom vcf-sort input.vcf > sorted.vcf
```

### Combined workflow

```bash
docker run -v $(pwd):/data biopsyk/custom bash -c '\
  bcftools view -i "QUAL>30" input.vcf.gz -Ov -o filtered.vcf && \
  vcf-sort filtered.vcf > sorted.vcf && \
  vcftools --vcf sorted.vcf --freq --out freqs'
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/custom:1.0.0

singularity exec custom_1.0.0.sif bcftools --help
singularity exec custom_1.0.0.sif vcftools --help
singularity exec custom_1.0.0.sif perl -v
```

## When to use standalone instead

| Need | Image |
|------|--------|
| Only BCFtools / tabix / bgzip | `biopsyk/bcftools` |
| Only a Perl interpreter | `biopsyk/perl` |
| Only VCFtools | `biopsyk/vcftools` |
| All three in one container | `biopsyk/custom` |

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- HTSlib `tabix`/`bgzip` come from the BCFtools build (1.22.1), not from the Ubuntu `tabix` package used in the standalone VCFtools image

## Building the Image

From this directory:

```bash
docker build -t biopsyk/custom:latest .
```

## References

- [BCFtools documentation](https://samtools.github.io/bcftools/bcftools.html)
- [Perl documentation](https://perldoc.perl.org/)
- [VCFtools website](https://vcftools.github.io/)
