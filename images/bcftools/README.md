# BCFtools Docker Image

This Docker image provides BCFtools along with HTSlib utilities (tabix, bgzip, htsfile) for manipulation and analysis of VCF/BCF files.

## Description

BCFtools is a set of utilities for variant calling and manipulating VCF (Variant Call Format) and BCF (binary VCF) files. It is part of the SAMtools project and provides tools for filtering, subsetting, concatenating, annotating, and converting variant files.

HTSlib provides low-level APIs for working with high-throughput sequencing data formats including SAM, BAM, CRAM, VCF, and BCF. This image includes command-line utilities from HTSlib:
- **tabix**: Indexes and queries TAB-delimited genome position files
- **bgzip**: Block-compressed GZIP format for genomic data
- **htsfile**: Identifies high-throughput sequencing data file formats

## Version Information

- Internal image version: 1.0.0
- BCFtools version: 1.22
- HTSlib version: 1.22.1

## Key Features

### BCFtools Capabilities
- **Filtering**: Filter variants by quality, depth, allele frequency, etc.
- **Subsetting**: Extract specific samples or regions
- **Merging**: Combine multiple VCF/BCF files
- **Annotation**: Add or modify variant annotations
- **Statistics**: Generate summary statistics from VCF/BCF files
- **Normalization**: Left-align and normalize indels
- **Conversion**: Convert between VCF and BCF formats
- **Consensus**: Create consensus sequences from VCF files

### HTSlib Tools
- **tabix**: Fast random access to compressed genomic data files
- **bgzip**: Compress/decompress files while maintaining seekability
- **htsfile**: File format detection and validation

## Usage

### Basic BCFtools Usage

```bash
# View VCF/BCF file
docker run -v $(pwd):/data biopsyk/bcftools bcftools view input.vcf.gz

# Filter variants by quality
docker run -v $(pwd):/data biopsyk/bcftools bcftools view -i 'QUAL>30' input.vcf.gz -o filtered.vcf.gz

# Extract specific samples
docker run -v $(pwd):/data biopsyk/bcftools bcftools view -s sample1,sample2 input.vcf.gz -o subset.vcf.gz

# Extract specific region
docker run -v $(pwd):/data biopsyk/bcftools bcftools view -r chr1:1000000-2000000 input.vcf.gz

# Get statistics
docker run -v $(pwd):/data biopsyk/bcftools bcftools stats input.vcf.gz > stats.txt

# Normalize variants (left-align and split multiallelic sites)
docker run -v $(pwd):/data biopsyk/bcftools bcftools norm -m-both -f reference.fa input.vcf.gz -o normalized.vcf.gz
```

### HTSlib Tools Usage

```bash
# Compress VCF with bgzip
docker run -v $(pwd):/data biopsyk/bcftools bgzip input.vcf

# Index VCF file with tabix
docker run -v $(pwd):/data biopsyk/bcftools tabix -p vcf input.vcf.gz

# Query specific region from indexed VCF
docker run -v $(pwd):/data biopsyk/bcftools tabix input.vcf.gz chr1:1000000-2000000

# Check file format
docker run -v $(pwd):/data biopsyk/bcftools htsfile input.vcf.gz
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/bcftools:1.0.0

# Run bcftools
singularity exec bcftools_1.0.0.sif bcftools view input.vcf.gz

# Run tabix
singularity exec bcftools_1.0.0.sif tabix -p vcf input.vcf.gz

# Run bgzip
singularity exec bcftools_1.0.0.sif bgzip input.vcf
```

## Common Analysis Workflows

### 1. Quality Control and Filtering

```bash
# Get basic statistics
docker run -v $(pwd):/data biopsyk/bcftools bcftools stats input.vcf.gz > stats.txt

# Filter by quality, depth, and allele frequency
docker run -v $(pwd):/data biopsyk/bcftools bcftools view \
  -i 'QUAL>30 && DP>10 && AF>0.01' \
  input.vcf.gz -Oz -o filtered.vcf.gz

# Index the filtered file
docker run -v $(pwd):/data biopsyk/bcftools tabix -p vcf filtered.vcf.gz
```

### 2. Merging VCF Files

```bash
# Merge multiple VCF files (same samples, different regions)
docker run -v $(pwd):/data biopsyk/bcftools bcftools concat \
  file1.vcf.gz file2.vcf.gz file3.vcf.gz \
  -Oz -o merged.vcf.gz

# Merge multiple VCF files (different samples, same regions)
docker run -v $(pwd):/data biopsyk/bcftools bcftools merge \
  sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
  -Oz -o merged.vcf.gz
```

### 3. Annotation and Modification

```bash
# Annotate variants with information from another VCF
docker run -v $(pwd):/data biopsyk/bcftools bcftools annotate \
  -a annotations.vcf.gz -c INFO/AF,INFO/AC \
  input.vcf.gz -Oz -o annotated.vcf.gz

# Add or modify tags
docker run -v $(pwd):/data biopsyk/bcftools bcftools annotate \
  --set-id '%CHROM\_%POS\_%REF\_%ALT' \
  input.vcf.gz -Oz -o with_ids.vcf.gz
```

### 4. Sample Manipulation

```bash
# Extract specific samples
docker run -v $(pwd):/data biopsyk/bcftools bcftools view \
  -s sample1,sample2,sample3 \
  input.vcf.gz -Oz -o subset.vcf.gz

# Remove specific samples
docker run -v $(pwd):/data biopsyk/bcftools bcftools view \
  -s ^sample_to_remove \
  input.vcf.gz -Oz -o without_sample.vcf.gz

# Rename samples
docker run -v $(pwd):/data biopsyk/bcftools bcftools reheader \
  -s sample_names.txt \
  input.vcf.gz -o renamed.vcf.gz
```

### 5. Variant Normalization

```bash
# Left-align and normalize indels, split multiallelic sites
docker run -v $(pwd):/data biopsyk/bcftools bcftools norm \
  -f reference.fa -m-both \
  input.vcf.gz -Oz -o normalized.vcf.gz

# Remove duplicates
docker run -v $(pwd):/data biopsyk/bcftools bcftools norm \
  -d all \
  input.vcf.gz -Oz -o dedup.vcf.gz
```

### 6. Format Conversion

```bash
# Convert VCF to BCF (binary, more efficient)
docker run -v $(pwd):/data biopsyk/bcftools bcftools view \
  input.vcf.gz -Ob -o output.bcf

# Convert BCF to VCF
docker run -v $(pwd):/data biopsyk/bcftools bcftools view \
  input.bcf -Oz -o output.vcf.gz

# Index BCF file
docker run -v $(pwd):/data biopsyk/bcftools bcftools index output.bcf
```

### 7. Consensus Sequence Generation

```bash
# Create consensus FASTA from reference and variants
docker run -v $(pwd):/data biopsyk/bcftools bcftools consensus \
  -f reference.fa \
  variants.vcf.gz > consensus.fa

# Create consensus for specific sample
docker run -v $(pwd):/data biopsyk/bcftools bcftools consensus \
  -f reference.fa -s sample1 \
  variants.vcf.gz > sample1_consensus.fa
```

## Input/Output Format Options

BCFtools supports multiple output formats via the `-O` flag:
- `-Ov`: Uncompressed VCF
- `-Oz`: Compressed VCF (gzip)
- `-Ob`: Compressed BCF (binary VCF)
- `-Ou`: Uncompressed BCF (for piping between bcftools commands)

Example of piping commands:
```bash
docker run -v $(pwd):/data biopsyk/bcftools bash -c \
  "bcftools view -i 'QUAL>30' input.vcf.gz -Ou | \
   bcftools norm -m-both -f ref.fa -Oz -o output.vcf.gz"
```

## Performance Tips

- Use BCF format (`-Ob`) instead of VCF for intermediate files - it's faster to read/write
- Use `-Ou` (uncompressed BCF) when piping between bcftools commands
- Always index your VCF/BCF files with `tabix` or `bcftools index` for random access
- For large files, specify regions (`-r` or `-R`) to avoid processing entire files
- Use multi-threading with `--threads` option when available

## File Requirements

### VCF/BCF Files
- Should be compressed with bgzip for indexing
- Must be indexed with tabix for random access queries
- Headers should contain proper contig definitions

### Reference FASTA
- Required for normalization and consensus calling
- Must be indexed with `samtools faidx`

## Common Filter Expressions

```bash
# Filter by quality
-i 'QUAL>30'

# Filter by depth
-i 'DP>10'

# Filter by allele frequency
-i 'AF>0.01'

# Combine multiple filters
-i 'QUAL>30 && DP>10 && AF>0.01'

# Filter by variant type
-i 'TYPE="snp"'           # SNPs only
-i 'TYPE="indel"'         # Indels only

# Filter biallelic sites only
-i 'N_ALT=1'

# Filter by missing data
-i 'F_MISSING<0.1'        # Less than 10% missing
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- VCF files should be compressed with bgzip and indexed with tabix for optimal performance
- The image includes all BCFtools plugins

## Building the Image

From this directory:

```bash
docker build -t biopsyk/bcftools:latest .
```

## References

- [BCFtools Documentation](https://samtools.github.io/bcftools/bcftools.html)
- [BCFtools GitHub Repository](https://github.com/samtools/bcftools)
- [HTSlib Documentation](http://www.htslib.org/)
- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- Danecek P, et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience, 10(2):giab008.
