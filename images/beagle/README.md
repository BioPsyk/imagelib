# Beagle Docker Image

This Docker image provides Beagle 5.5, a software package for phasing genotypes and imputing ungenotyped markers.

## Description

Beagle is a state-of-the-art software package for analyzing large-scale genetic data. It performs genotype phasing, genotype imputation, and genetic association analysis. Beagle is designed to be fast and memory-efficient, capable of handling millions of genetic variants across thousands of samples.

## Version Information

- Internal image version: 1.0.0
- Beagle version: 5.5 (27Feb25.75f)
- Built from: https://faculty.washington.edu/browning/beagle/

## Key Features

- **Genotype Phasing**: Estimates haplotype phase from unphased genotype data
- **Genotype Imputation**: Imputes ungenotyped markers using a reference panel
- **High Performance**: Optimized algorithms for speed and memory efficiency
- **Large-Scale Analysis**: Handles millions of markers and thousands of samples
- **Reference Panel Support**: Works with various reference panels (1000 Genomes, TOPMed, HRC, etc.)
- **VCF Format**: Reads and writes standard VCF and BCF files
- **Identity by Descent (IBD)**: Detects IBD segments between samples
- **Genetic Association**: Performs single-marker and multi-marker association tests

## Usage

### Basic Commands

```bash
# Show help and version
docker run biopsyk/beagle

# Run Beagle with parameters
docker run -v $(pwd):/data biopsyk/beagle beagle gt=input.vcf.gz out=output

# With custom memory allocation (default is 4GB)
docker run -v $(pwd):/data biopsyk/beagle bash -c "java -Xmx16g -jar /opt/beagle/beagle.27Feb25.75f.jar gt=input.vcf.gz out=output"
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/beagle:1.0.0

# Run Beagle
singularity exec beagle_1.0.0.sif beagle gt=input.vcf.gz out=output

# With custom memory allocation
singularity exec beagle_1.0.0.sif bash -c "java -Xmx16g -jar /opt/beagle/beagle.27Feb25.75f.jar gt=input.vcf.gz out=output"
```

## Common Workflows

### 1. Genotype Phasing

Phase genotypes without imputation:

```bash
# Basic phasing
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=unphased.vcf.gz out=phased

# Phasing with a genetic map
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=unphased.vcf.gz map=plink.chr22.map out=phased
```

Output files:
- `phased.vcf.gz`: Phased genotypes in VCF format
- `phased.log`: Log file with analysis details

### 2. Genotype Imputation

Impute ungenotyped markers using a reference panel:

```bash
# Basic imputation
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=target.vcf.gz ref=reference.vcf.gz out=imputed

# Imputation with genetic map
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=target.vcf.gz ref=reference.vcf.gz map=plink.map out=imputed

# Imputation with quality filtering
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=target.vcf.gz ref=reference.vcf.gz out=imputed impute=true gp=true
```

Parameters:
- `gt=`: Target genotypes to be phased and imputed
- `ref=`: Reference panel for imputation
- `map=`: Genetic map file
- `out=`: Output file prefix
- `impute=true`: Enable imputation (default)
- `gp=true`: Output genotype probabilities

### 3. Phasing and Imputation with a Reference Panel

```bash
# Combined phasing and imputation
docker run -v $(pwd):/data biopsyk/beagle \
  beagle \
  gt=target.vcf.gz \
  ref=reference.vcf.gz \
  map=genetic_map.txt \
  out=phased_imputed \
  nthreads=8
```

### 4. IBD Detection

Detect identity-by-descent (IBD) segments:

```bash
# Basic IBD detection
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=phased.vcf.gz out=ibd ibd=true

# IBD with custom LOD threshold
docker run -v $(pwd):/data biopsyk/beagle \
  beagle gt=phased.vcf.gz out=ibd ibd=true ibdlod=3.0
```

Output files:
- `ibd.ibd.gz`: IBD segments
- `ibd.hbd.gz`: Homozygosity-by-descent segments

### 5. Imputation with 1000 Genomes Reference

Example using 1000 Genomes Phase 3 as reference:

```bash
# Download 1000G reference panel (example for chromosome 22)
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz

# Download genetic map
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.chr22.GRCh37.map

# Run imputation
docker run -v $(pwd):/data biopsyk/beagle \
  beagle \
  gt=target.chr22.vcf.gz \
  ref=chr22.1kg.phase3.v5a.vcf.gz \
  map=plink.chr22.GRCh37.map \
  out=imputed.chr22 \
  nthreads=4
```

### 6. Quality Control and Filtering

```bash
# Filter by INFO score after imputation
docker run -v $(pwd):/data biopsyk/beagle \
  beagle \
  gt=target.vcf.gz \
  ref=reference.vcf.gz \
  out=imputed \
  gprobs=true

# Then filter using bcftools (if available)
bcftools view -i 'INFO/DR2>0.8' imputed.vcf.gz -Oz -o imputed.filtered.vcf.gz
```

## Beagle Parameters Reference

### Input/Output Parameters

```
gt=<file>           Input VCF file with target genotypes
ref=<file>          Reference panel VCF file (for imputation)
out=<prefix>        Output file prefix (default: beagle)
map=<file>          Genetic map file (PLINK format)
excludesamples=<file>  File with samples to exclude
excludemarkers=<file>  File with markers to exclude
```

### Phasing and Imputation Parameters

```
impute=<true/false>     Impute ungenotyped markers (default: true)
gp=<true/false>         Output genotype probabilities (default: false)
ap=<true/false>         Output allele probabilities (default: false)
niterations=<int>       Number of phasing iterations (default: 5)
phase-states=<int>      Number of phasing states (default: 280)
impute-states=<int>     Number of imputation states (default: 1600)
burnin=<int>            Number of burnin iterations (default: 3)
```

### IBD Detection Parameters

```
ibd=<true/false>        Detect IBD segments (default: false)
ibdlod=<float>          LOD score threshold for reporting IBD (default: 3.0)
ibdcm=<float>           Minimum cM length for reporting IBD (default: 2.0)
```

### Performance Parameters

```
nthreads=<int>          Number of threads (default: all available)
window=<float>          Window size in cM for analysis (default: 40.0)
overlap=<float>         Window overlap in cM (default: 4.0)
```

### Quality Control Parameters

```
lowmem=<true/false>     Use low memory mode (default: false)
gprobs=<true/false>     Output genotype probabilities (default: false)
ne=<int>                Effective population size (default: 1000000)
err=<float>             Allele mismatch rate (default: 0.0001)
```

## Input File Format

### Target Genotypes (VCF Format)

Beagle accepts standard VCF or BCF files:

```
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID         REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
22      16050075 rs587697622  A    G    .     PASS    .     GT      0/1      0/0
22      16050115 rs587755077  G    A    .     PASS    .     GT      1/1      0/1
```

Requirements:
- Must contain GT (genotype) field
- Can be unphased (0/1) or phased (0|1)
- Can be compressed (.vcf.gz) or uncompressed (.vcf)
- BCF format (.bcf) is also supported

### Genetic Map Format

PLINK format genetic map:

```
chr position rate map
22 16050075 3.79978 0.000000
22 16050115 3.79978 0.000152
22 16050213 3.79978 0.000524
```

Columns:
1. Chromosome
2. Physical position (bp)
3. Recombination rate (cM/Mb)
4. Genetic position (cM)

## Output Files

### Phased/Imputed Genotypes

Main output file: `<prefix>.vcf.gz`

Contains:
- **GT**: Phased genotypes (0|1, 1|0, etc.)
- **DS**: Dosage (allele count: 0, 1, or 2)
- **GP**: Genotype probabilities (if `gp=true`)
- **DR2**: Dosage r² (imputation quality metric)

### Log File

File: `<prefix>.log`

Contains:
- Command line used
- Data summary statistics
- Runtime and memory usage
- Number of markers and samples processed
- Imputation quality summaries

### IBD Output Files

**IBD segments**: `<prefix>.ibd.gz`

Format:
```
sample1  sample2  chrom  start_pos  end_pos  LOD  length_cM
```

**HBD segments**: `<prefix>.hbd.gz`

Format:
```
sample  chrom  start_pos  end_pos  LOD  length_cM
```

## Interpreting Results

### Imputation Quality (DR2)

The DR2 (dosage r²) metric measures imputation accuracy:

- **DR2 > 0.8**: High quality imputation, suitable for association analysis
- **DR2 0.5-0.8**: Moderate quality, use with caution
- **DR2 < 0.5**: Low quality, consider excluding from analysis

### Genotype Probabilities (GP)

When `gp=true`, output includes three probabilities for each genotype:
- P(0/0): Probability of homozygous reference
- P(0/1): Probability of heterozygous
- P(1/1): Probability of homozygous alternate

Example: `GP=0.98,0.02,0.00` indicates 98% probability of 0/0 genotype.

### IBD Segments

- **LOD score**: Log-odds score for IBD vs. non-IBD
  - LOD > 3.0: Strong evidence for IBD
  - LOD > 5.0: Very strong evidence
- **Length in cM**: Genetic length of IBD segment
  - Recent common ancestry: longer segments
  - Distant common ancestry: shorter segments

## Best Practices

### 1. Input Data Preparation

- **Quality Control**: Filter variants and samples before phasing/imputation
- **Format**: Use VCF format with GT field
- **Compression**: Compress large files with bgzip
- **Chromosome Splitting**: Process chromosomes separately for efficiency

### 2. Reference Panel Selection

- **Match Ancestry**: Use reference panels matching your study population
- **Size**: Larger panels generally improve imputation accuracy
- **Common Panels**:
  - 1000 Genomes Phase 3 (2,504 samples, all ancestries)
  - TOPMed (>130,000 samples, primarily European/African)
  - HRC (32,470 samples, primarily European)

### 3. Memory Management

- **Default**: 4GB allocated by wrapper script
- **Large Datasets**: Increase with `-Xmx` flag
  - 1000 samples: 8GB recommended
  - 10,000 samples: 16-32GB recommended
  - 50,000+ samples: 64GB+ may be needed

### 4. Performance Optimization

- **Threading**: Use `nthreads=` parameter to utilize multiple cores
- **Window Size**: Adjust `window=` for balance between speed and accuracy
- **Low Memory Mode**: Use `lowmem=true` for very large datasets

### 5. Quality Filtering

- **Pre-imputation**: Remove low-quality variants from target data
- **Post-imputation**: Filter by DR2 score (typically DR2 > 0.3 or 0.8)
- **MAF**: Consider minor allele frequency filters

## Common Issues and Solutions

### Issue: "OutOfMemoryError: Java heap space"
**Solution**: Increase memory allocation using `-Xmx` flag:
```bash
docker run -v $(pwd):/data biopsyk/beagle bash -c "java -Xmx16g -jar /opt/beagle/beagle.27Feb25.75f.jar gt=input.vcf.gz out=output"
```

### Issue: "Inconsistent markers in target and reference"
**Solution**: Ensure chromosome names match between files. Use bcftools to rename:
```bash
bcftools annotate --rename-chrs chr_name_conv.txt input.vcf.gz -Oz -o output.vcf.gz
```

### Issue: Slow performance
**Solution**:
1. Use `nthreads=` to enable parallel processing
2. Process chromosomes separately
3. Use `lowmem=true` if memory constrained

### Issue: Low imputation accuracy
**Solution**:
1. Use larger or better-matched reference panel
2. Ensure genetic map is correct for your genome build
3. Improve quality control of input data

## Performance Notes

- **Speed**: Processes ~1M markers in 10-30 minutes (depends on sample size)
- **Memory**: Scales with sample size and marker density
- **Threading**: Scales well with multiple cores
- **Disk I/O**: Bottleneck for very large files, use compression

## Example Complete Workflow

```bash
# Step 1: Quality control of input VCF
bcftools view -i 'F_MISSING<0.1 && MAF>0.01' raw.vcf.gz -Oz -o qc.vcf.gz
bcftools index qc.vcf.gz

# Step 2: Download reference panel and map
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.chr22.GRCh37.map

# Step 3: Run Beagle imputation
docker run -v $(pwd):/data biopsyk/beagle \
  beagle \
  gt=qc.vcf.gz \
  ref=chr22.1kg.phase3.v5a.vcf.gz \
  map=plink.chr22.GRCh37.map \
  out=imputed.chr22 \
  nthreads=8 \
  gp=true

# Step 4: Filter by imputation quality
bcftools view -i 'DR2>0.8' imputed.chr22.vcf.gz -Oz -o imputed.chr22.filtered.vcf.gz
bcftools index imputed.chr22.filtered.vcf.gz

# Results:
# - imputed.chr22.filtered.vcf.gz: High-quality imputed genotypes
# - imputed.chr22.log: Analysis log with quality metrics
```

## Building the Image

From this directory:

```bash
docker build -t biopsyk/beagle:latest .
```

## References

### Primary Citations

- Browning, B.L., Tian, X., Zhou, Y., and Browning, S.R. (2021). Fast two-stage phasing of large-scale sequence data. *American Journal of Human Genetics*, 108(10):1880-1890. [doi:10.1016/j.ajhg.2021.08.005](https://doi.org/10.1016/j.ajhg.2021.08.005)

- Browning, B.L., Zhou, Y., and Browning, S.R. (2018). A one-penny imputed genome from next-generation reference panels. *American Journal of Human Genetics*, 103(3):338-348. [doi:10.1016/j.ajhg.2018.07.015](https://doi.org/10.1016/j.ajhg.2018.07.015)

### Additional Resources

- [Beagle 5.5 Documentation](https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf)
- [Beagle Official Website](https://faculty.washington.edu/browning/beagle/beagle.html)
- [1000 Genomes Reference Panels](http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/)
- [Genetic Maps](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

## Notes

- This image runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- Default memory allocation is 4GB; increase for large datasets
- Java 21 is used for optimal performance

## License

Beagle is distributed under the GNU General Public License version 3 (GPLv3). See the [LICENSE](https://www.gnu.org/licenses/gpl-3.0.html) for details.

## Authors

Developed by:
- Brian L. Browning (University of Washington)
- Sharon R. Browning (University of Washington)
