# Perl Docker Image

This Docker image provides a full Perl 5 interpreter for running standalone Perl scripts.

## Description

Perl is a general-purpose scripting language commonly used in bioinformatics for text processing, format conversion, and glue code around other tools.

## Version Information

- Internal image version: 1.0.0
- Perl version: 5.34 (Ubuntu 22.04)

The installed version can be verified by running `perl -v` or by reading `/etc/perl-version` inside the container.

## Usage

### Basic Usage

```bash
# Show Perl version
docker run -v $(pwd):/data biopsyk/perl perl -v

# Run a script mounted from the host
docker run -v $(pwd):/data biopsyk/perl perl myscript.pl

# One-liner
docker run -v $(pwd):/data biopsyk/perl perl -ne 'print if /rs[0-9]+/' input.txt
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/perl:1.0.0

# Run a script
singularity exec perl_1.0.0.sif perl myscript.pl
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- This image contains the Perl interpreter only. VCFtools Perl helpers (`vcf-sort`, `Vcf.pm`, etc.) live in the `vcftools` and `custom` images.

## Building the Image

From this directory:

```bash
docker build -t biopsyk/perl:latest .
```

## References

- [Perl documentation](https://perldoc.perl.org/)
- [Perl downloads](https://www.perl.org/get.html)
