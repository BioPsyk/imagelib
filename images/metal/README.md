# METAL Docker Image

This Docker image provides METAL, a tool for meta-analysis of genome-wide association studies (GWAS).

## Description

METAL (Meta-Analysis Helper) is a fast and efficient tool for meta-analysis of genomewide association scans. It implements various meta-analysis approaches including sample-size based, inverse-variance based, and custom meta-analysis schemes. METAL can handle millions of genetic variants across multiple studies efficiently.

## Version Information

- Internal image version: 1.0.0
- METAL version: Latest from GitHub (master branch)
- Built from: https://github.com/statgen/METAL

## Key Features

- **Multiple Meta-Analysis Schemes**: Sample-size weighted, inverse-variance weighted, and custom approaches
- **Efficient Processing**: Handles millions of variants across multiple studies
- **Genomic Control**: Built-in genomic control correction
- **Flexible Input**: Supports various GWAS summary statistic formats
- **Heterogeneity Testing**: Cochran's Q statistic for heterogeneity assessment
- **Sample Overlap Correction**: Accounts for overlapping samples across studies
- **Position Tracking**: Chromosome and position validation across studies
- **Allele Frequency Analysis**: Meta-analysis of allele frequencies

## Usage

### Basic Commands

```bash
# Show help
docker run -v $(pwd):/data biopsyk/metal metal --help

# Run METAL with a script file
docker run -v $(pwd):/data biopsyk/metal metal < metal_script.txt
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/metal:1.0.0

# Run METAL
singularity exec metal_1.0.0.sif metal --help
```

## Common Meta-Analysis Workflows

### 1. Sample-Size Weighted Meta-Analysis

This is the simplest approach, weighting studies by sample size.

Create a METAL script file (`metal_samplesize.txt`):

```
# Sample size weighted meta-analysis
SCHEME SAMPLESIZE

# === DESCRIBE AND PROCESS FIRST FILE ===
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
PVALUE   P
WEIGHT   N

PROCESS study1_results.txt

# === DESCRIBE AND PROCESS SECOND FILE ===
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
PVALUE   P
WEIGHT   N

PROCESS study2_results.txt

# === DESCRIBE AND PROCESS THIRD FILE ===
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
PVALUE   P
WEIGHT   N

PROCESS study3_results.txt

# === OUTPUT OPTIONS ===
OUTFILE meta_analysis_results .tbl

# === RUN ANALYSIS ===
ANALYZE HETEROGENEITY

QUIT
```

Run the analysis:

```bash
docker run -v $(pwd):/data biopsyk/metal metal < metal_samplesize.txt
```

### 2. Inverse-Variance Weighted Meta-Analysis

More powerful approach using effect sizes and standard errors.

Create a METAL script file (`metal_stderr.txt`):

```
# Inverse variance weighted meta-analysis
SCHEME STDERR

# === GENOMIC CONTROL ===
GENOMICCONTROL ON

# === DESCRIBE AND PROCESS FILES ===
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
STDERR   SE
PVALUE   P

PROCESS study1_results.txt

MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
STDERR   SE
PVALUE   P

PROCESS study2_results.txt

MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
STDERR   SE
PVALUE   P

PROCESS study3_results.txt

# === OUTPUT OPTIONS ===
OUTFILE meta_analysis_stderr .tbl

# === RUN ANALYSIS ===
ANALYZE HETEROGENEITY

QUIT
```

Run the analysis:

```bash
docker run -v $(pwd):/data biopsyk/metal metal < metal_stderr.txt
```

### 3. Meta-Analysis with Sample Overlap Correction

When studies have overlapping samples:

```
# Sample size weighted with overlap correction
SCHEME SAMPLESIZE
OVERLAP ON
ZCUTOFF 1.0

# === DESCRIBE AND PROCESS FILES ===
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA
PVALUE   P
WEIGHT   N

PROCESS study1_results.txt
PROCESS study2_results.txt
PROCESS study3_results.txt

# === OUTPUT OPTIONS ===
OUTFILE meta_overlap_corrected .tbl

# === RUN ANALYSIS ===
ANALYZE HETEROGENEITY

QUIT
```

### 4. Position Tracking for Manhattan Plots

Track chromosome and position information:

```
# Enable position tracking
SCHEME STDERR
TRACKPOSITIONS ON

MARKER     SNP
ALLELE     A1 A2
EFFECT     BETA
STDERR     SE
PVALUE     P
CHROMOSOME CHR
POSITION   BP

PROCESS study1_results.txt
PROCESS study2_results.txt

OUTFILE meta_with_positions .tbl

ANALYZE HETEROGENEITY

QUIT
```

### 5. Custom Analysis with Multiple Trait Meta-Analysis

For correlated traits or multi-phenotype analysis:

```
SCHEME STDERR

# Analyze trait 1
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA_trait1
STDERR   SE_trait1

PROCESS trait1_study1.txt
PROCESS trait1_study2.txt

OUTFILE trait1_meta .tbl
ANALYZE

CLEAR

# Analyze trait 2
MARKER   SNP
ALLELE   A1 A2
EFFECT   BETA_trait2
STDERR   SE_trait2

PROCESS trait2_study1.txt
PROCESS trait2_study2.txt

OUTFILE trait2_meta .tbl
ANALYZE

QUIT
```

## METAL Script Commands Reference

### Analysis Schemes

```
SCHEME SAMPLESIZE    # Sample-size weighted
SCHEME STDERR        # Inverse-variance weighted (default)
```

### Input File Description

```
MARKER    column_name   # SNP identifier column
ALLELE    col1 col2     # Effect and non-effect allele columns
EFFECT    column_name   # Effect size (BETA or OR)
STDERR    column_name   # Standard error (required for STDERR scheme)
PVALUE    column_name   # P-value column
WEIGHT    column_name   # Sample size (required for SAMPLESIZE scheme)
CHROMOSOME column_name  # Chromosome column (for TRACKPOSITIONS)
POSITION   column_name  # Position column (for TRACKPOSITIONS)
FREQ      column_name   # Allele frequency column (optional)
```

### Processing Options

```
PROCESS filename          # Process input file
SEPARATOR WHITESPACE      # Tab/space delimited (default)
SEPARATOR COMMA           # Comma delimited
COLUMNCOUNTING STRICT     # Enforce column counting
GENOMICCONTROL ON         # Apply genomic control correction
OVERLAP ON                # Correct for sample overlap (SAMPLESIZE only)
TRACKPOSITIONS ON         # Track chromosome and position
CUSTOMVARIABLE variable   # Define custom variable
```

### Output Options

```
OUTFILE prefix suffix     # Output file naming
MAXWARNINGS n             # Maximum warnings to display
MINWEIGHT n               # Minimum weight threshold
EFFECT_PRINT_PRECISION n  # Digits for effect column (default: 4)
STDERR_PRINT_PRECISION n  # Digits for stderr column (default: 4)
```

### Analysis Commands

```
ANALYZE                   # Run meta-analysis
ANALYZE HETEROGENEITY     # Include heterogeneity statistics
CLEAR                     # Clear current analysis
QUIT                      # Exit METAL
```

## Input File Format

METAL accepts tab or space-delimited text files with headers. Example:

```
SNP          A1  A2  BETA      SE       P         N
rs1234567    A   G   0.0245    0.0123   0.0456    5000
rs2345678    C   T   -0.0312   0.0145   0.0321    5000
rs3456789    G   A   0.0189    0.0134   0.1234    5000
```

### Required Columns

For **STDERR** scheme:
- SNP identifier (MARKER)
- Alleles (ALLELE)
- Effect size (EFFECT)
- Standard error (STDERR)

For **SAMPLESIZE** scheme:
- SNP identifier (MARKER)
- Alleles (ALLELE)
- Effect size (EFFECT)
- P-value (PVALUE)
- Sample size (WEIGHT)

### Optional Columns
- Chromosome (CHR)
- Position (BP/POS)
- Allele frequency (FREQ/MAF)
- Effect allele frequency (EAF)

## Output Files

METAL generates several output files:

### Main Results File (.tbl)

Columns in output:
- **MarkerName**: SNP identifier
- **Allele1**: Effect allele
- **Allele2**: Non-effect allele
- **Freq1**: Meta-analyzed allele frequency
- **Effect**: Meta-analyzed effect size
- **StdErr**: Meta-analyzed standard error
- **P-value**: Meta-analysis P-value
- **Direction**: Direction of effects across studies (e.g., "+++", "+--")
- **HetISq**: I² heterogeneity statistic
- **HetChiSq**: Cochran's Q statistic
- **HetDf**: Degrees of freedom for heterogeneity test
- **HetPVal**: P-value for heterogeneity

With TRACKPOSITIONS:
- **Chromosome**: Chromosome
- **Position**: Base-pair position

### Additional Output Files

- `.info`: Information about each study processed
- `.log`: Log file with analysis details
- `.txt`: Additional results (if specified)

## Interpreting Results

### Effect Direction String

The Direction column shows effect direction across studies:
- `+++`: All three studies show positive effect
- `++-`: Two positive, one negative
- `+--`: One positive, two negative
- `?+-`: One study missing data

### Heterogeneity Statistics

- **HetISq (I²)**: Percentage of variation due to heterogeneity
  - 0-25%: Low heterogeneity
  - 25-50%: Moderate heterogeneity
  - 50-75%: Substantial heterogeneity
  - 75-100%: Considerable heterogeneity

- **HetPVal < 0.05**: Significant heterogeneity detected

### Genomic Control

When enabled, METAL calculates genomic inflation factor (λ):
- λ ≈ 1.0: No inflation
- λ > 1.05: May indicate population stratification
- Effect sizes and standard errors are corrected by √λ

## Best Practices

1. **Quality Control**:
   - Pre-filter studies for INFO > 0.3 or 0.4
   - Remove SNPs with very low MAF
   - Check for strand issues

2. **Consistent Effect Coding**:
   - Ensure effect alleles are consistent across studies
   - METAL will flip alleles automatically when possible
   - Check Direction column for unexpected flips

3. **Genomic Control**:
   - Apply to individual studies before meta-analysis
   - Or apply within METAL using GENOMICCONTROL ON
   - Don't apply twice (study-level + METAL)

4. **Heterogeneity**:
   - Always analyze heterogeneity (ANALYZE HETEROGENEITY)
   - Investigate SNPs with high I² or low HetPVal
   - Consider fixed vs random effects

5. **Sample Overlap**:
   - Use OVERLAP ON if studies share participants
   - Document overlap percentages
   - Consider sensitivity analyses

6. **Output Precision**:
   - Increase EFFECT_PRINT_PRECISION for downstream use
   - Default 4 decimals may lose precision for small effects

## Common Issues and Solutions

### Issue: "Alleles do not match"
**Solution**: Check that allele coding is consistent. METAL will attempt to flip alleles, but review warnings.

### Issue: "No overlapping markers"
**Solution**: Ensure SNP identifiers are consistent across studies (rsIDs or chr:pos format).

### Issue: High genomic inflation (λ > 1.1)
**Solution**: Check for population stratification. May need to apply stricter QC or adjust analysis.

### Issue: Very different effect sizes across studies
**Solution**: Check for heterogeneity. May indicate true biological differences or technical issues.

## Performance Notes

- METAL is very fast, processing millions of SNPs in minutes
- Memory usage is low (typically < 1GB)
- Can handle 10+ studies easily
- No practical limit on number of SNPs

## Example Run

```bash
# Create a simple meta-analysis script
cat > metal_script.txt << 'EOF'
SCHEME STDERR
GENOMICCONTROL ON

MARKER SNP
ALLELE A1 A2
EFFECT BETA
STDERR SE
PVALUE P

PROCESS study1.txt
PROCESS study2.txt
PROCESS study3.txt

OUTFILE myMetaAnalysis .tbl
ANALYZE HETEROGENEITY

QUIT
EOF

# Run METAL
docker run -v $(pwd):/data biopsyk/metal metal < metal_script.txt

# Results will be in myMetaAnalysis1.tbl
```

## Building the Image

From this directory:

```bash
docker build -t biopsyk/metal:latest .
```

## References

### Primary Citation

- Willer CJ, Li Y, Abecasis GR. (2010). METAL: fast and efficient meta-analysis of genomewide association scans. *Bioinformatics*, 26(17):2190-2191. [doi:10.1093/bioinformatics/btq340](https://doi.org/10.1093/bioinformatics/btq340)

### Additional Resources

- [METAL Documentation](https://genome.sph.umich.edu/wiki/METAL_Documentation)
- [METAL GitHub Repository](https://github.com/statgen/METAL)
- [METAL Official Website](http://csg.sph.umich.edu/abecasis/metal/)

## Notes

- This image runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- METAL reads commands from stdin, so use `< script.txt` or pipe commands
- The image includes the latest version from the GitHub master branch

## License

METAL is open source software. See the LICENSE file in the GitHub repository for details.

## Authors

Developed by:
- Cristen J. Willer (University of Michigan)
- Yun Li (University of North Carolina)
- Gonçalo R. Abecasis (University of Michigan)
