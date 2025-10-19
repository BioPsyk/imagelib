# LDSC Docker Image

This Docker image provides LDSC (LD Score Regression), a command-line tool for estimating heritability and genetic correlation from GWAS summary statistics.

## Description

LDSC (LD Score Regression) implements a method that distinguishes confounding biases from true polygenic signals in genome-wide association studies (GWAS). By regressing GWAS test statistics on linkage disequilibrium (LD) scores, LDSC estimates the proportion of phenotypic variance explained by genetic factors (heritability) and assesses the shared genetic architecture between traits (genetic correlation).

This image is based on the Python 3 implementation from [abdenlab/ldsc-python3](https://github.com/abdenlab/ldsc-python3), which modernizes the original LDSC tool for compatibility with current Python versions.

## Version Information

- Internal image version: 1.0.0
- LDSC version: 2.0.0 (Python 3 compatible)
- Python version: 3.10

## Key Features

- **Heritability Estimation**: Estimate SNP-based heritability from GWAS summary statistics
- **Genetic Correlation**: Calculate genetic correlations between traits
- **Partitioned Heritability**: Partition heritability by functional annotations
- **LD Score Computation**: Generate LD scores from genotype data
- **Regression Intercept**: Estimate and account for confounding biases
- **Summary Statistics Munging**: Standardize GWAS summary statistics for analysis

## Usage

### Basic Commands

```bash
# Show help for main LDSC tool
docker run -v $(pwd):/data biopsyk/ldsc ldsc -h

# Show help for summary statistics munging
docker run -v $(pwd):/data biopsyk/ldsc munge_sumstats -h

# Show help for annotation creation
docker run -v $(pwd):/data biopsyk/ldsc make_annot -h
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/ldsc:1.0.0

# Run LDSC
singularity exec ldsc_1.0.0.sif ldsc -h
```

## Common Analysis Workflows

### 1. Preparing Summary Statistics

Before running LDSC, you need to format (munge) your GWAS summary statistics:

```bash
docker run -v $(pwd):/data biopsyk/ldsc munge_sumstats \
  --sumstats mydata.sumstats.gz \
  --merge-alleles w_hm3.snplist \
  --out mydata \
  --a1-inc
```

### 2. Estimating Heritability

Estimate SNP-based heritability from GWAS summary statistics:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 mydata.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out mydata_h2
```

### 3. Calculating Genetic Correlation

Estimate genetic correlation between two traits:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --rg trait1.sumstats.gz,trait2.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out trait1_trait2_rg
```

### 4. Partitioned Heritability

Partition heritability by functional annotations:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 mydata.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out mydata_partitioned_h2
```

### 5. Computing LD Scores

Generate LD scores from your own genotype data:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --bfile mydata \
  --l2 \
  --ld-wind-cm 1 \
  --out mydata_ldscores
```

### 6. Cell Type-Specific Analysis

Test for enrichment in cell type-specific annotations:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2-cts mydata.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --ref-ld-chr-cts cell_type_data.ldcts \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --out mydata_celltype
```

## Stratified LD Score Regression (S-LDSC)

**Note**: S-LDSC functionality is built directly into LDSC and does not require a separate tool. The partitioned heritability and cell-type specific analyses shown above are implementations of S-LDSC.

### What is S-LDSC?

Stratified LD Score Regression (S-LDSC) extends basic LDSC to partition heritability across functional categories or annotations. This allows you to:
- Identify which genomic annotations are enriched for heritability
- Quantify the contribution of different functional categories to trait heritability
- Test hypotheses about the biological mechanisms underlying GWAS signals

### S-LDSC Workflow

#### Step 1: Create Custom Annotations

If you want to test your own genomic regions (e.g., enhancers, promoters, or tissue-specific regulatory elements):

```bash
# Create binary annotation file from BED file
docker run -v $(pwd):/data biopsyk/ldsc make_annot \
  --bed-file my_regions.bed \
  --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bim \
  --annot-file my_annotation.1.annot.gz
```

Repeat for all 22 chromosomes (replace `.1.` with chromosome numbers).

#### Step 2: Compute LD Scores for Custom Annotations

```bash
# Compute LD scores for your annotation
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --l2 \
  --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.1 \
  --ld-wind-cm 1 \
  --annot my_annotation.1.annot.gz \
  --out my_annotation.1
```

Repeat for all chromosomes.

#### Step 3: Run S-LDSC Analysis

Basic S-LDSC with single annotation category:

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 trait.sumstats.gz \
  --ref-ld-chr my_annotation. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out trait_sldsc_results
```

S-LDSC with baseline model (recommended):

```bash
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 trait.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline.,my_annotation. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out trait_sldsc_baseline_results
```

### S-LDSC with Baseline Annotations

The baseline model includes 53 functional categories (coding, UTR, promoter, enhancer, etc.):

```bash
# Download baseline annotations
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.2_ldscores.tgz
tar -xvzf 1000G_Phase3_baseline_v1.2_ldscores.tgz

# Run S-LDSC with baseline
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 trait.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out trait_baseline_sldsc \
  --print-coefficients
```

### S-LDSC for Tissue/Cell-Type Specific Analysis

Test enrichment in tissue-specific regulatory elements:

```bash
# Download tissue-specific annotations (e.g., GTEx)
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/Multi_tissue_gene_expr.ldcts

# Run cell-type specific S-LDSC
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2-cts trait.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --ref-ld-chr-cts Multi_tissue_gene_expr.ldcts \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --out trait_tissue_specific
```

### Interpreting S-LDSC Results

The output `.results` file contains:

- **Prop._SNPs**: Proportion of SNPs in the annotation
- **Prop._h2**: Proportion of heritability explained by the annotation
- **Prop._h2_std_error**: Standard error of proportion of heritability
- **Enrichment**: Fold-enrichment of heritability (Prop._h2 / Prop._SNPs)
- **Enrichment_std_error**: Standard error of enrichment
- **Enrichment_p**: P-value for enrichment test
- **Coefficient**: Regression coefficient (tau)
- **Coefficient_std_error**: Standard error of coefficient

**Key metrics**:
- **Enrichment > 1**: Annotation is enriched for heritability
- **Enrichment_p < 0.05**: Statistically significant enrichment
- **Prop._h2**: Shows biological importance (how much heritability is explained)

### S-LDSC Best Practices

1. **Always use baseline model**: Include `baseline_v1.2/baseline.` to control for known annotations
2. **Use --overlap-annot**: Accounts for overlap between annotations
3. **Use --print-coefficients**: Get regression coefficients (tau values)
4. **Multiple testing correction**: Apply Bonferroni or FDR correction for multiple annotations
5. **Check sample size**: Ensure GWAS has sufficient power (N > 5,000)
6. **Validate results**: Check if enrichment is robust to different LD reference panels

### Common S-LDSC Use Cases

**1. Test specific genomic regions**:
```bash
# E.g., test if brain enhancers are enriched for schizophrenia heritability
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 scz.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline.,brain_enhancers. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out scz_brain_enhancers
```

**2. Compare across tissues**:
```bash
# Use .ldcts file to test multiple tissues simultaneously
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2-cts trait.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --ref-ld-chr-cts GTEx_tissues.ldcts \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --out trait_gtex_tissues
```

**3. Custom functional categories**:
```bash
# Test your own ChIP-seq peaks, ATAC-seq peaks, etc.
docker run -v $(pwd):/data biopsyk/ldsc ldsc \
  --h2 trait.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline.,my_chipseq. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --out trait_chipseq_enrichment
```

### S-LDSC Resources

- [Baseline model annotations](https://data.broadinstitute.org/alkesgroup/LDSCORE/)
- [Cell-type specific annotations](https://alkesgroup.broadinstitute.org/LDSCORE/LDSC_SEG_ldscores/)
- [1000 Genomes reference panel](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz)
- [Weights for regression](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz)
- [Frequency files](https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz)

## Input File Formats

### GWAS Summary Statistics

LDSC requires summary statistics with the following columns:
- **SNP**: SNP identifier (rsID)
- **A1**: Effect allele
- **A2**: Non-effect allele
- **N**: Sample size
- **Z**: Z-score or **P**: P-value
- **OR** or **BETA**: Effect size

Optional columns:
- **INFO**: Imputation quality score
- **FRQ**: Allele frequency

### Reference Files

You'll need reference LD scores and weights. Download from:
- [European LD Scores](https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2)
- [East Asian LD Scores](https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2)
- [Partitioned LD Scores](https://data.broadinstitute.org/alkesgroup/LDSCORE/)

## Output Files

LDSC generates several output files:

### Heritability Output (.log file)
- Total SNP heritability estimate
- Standard error
- Lambda GC (genomic inflation factor)
- Mean chi-squared statistic
- LD Score regression intercept

### Genetic Correlation Output (.log file)
- Genetic correlation estimate (rg)
- Standard error
- P-value
- Heritability estimates for both traits

### Partitioned Heritability Output
- Heritability enrichment per annotation
- Proportion of heritability per annotation
- Coefficient and standard error

## Best Practices

1. **Quality Control**:
   - Remove SNPs with low imputation quality (INFO < 0.9)
   - Filter SNPs with extreme allele frequencies
   - Remove duplicate SNPs

2. **Sample Size**:
   - LDSC requires large GWAS (N > 5,000 recommended)
   - Ensure accurate sample size reporting

3. **LD Reference Panel**:
   - Use LD scores matching your study population ancestry
   - Use HapMap3 SNPs for stability

4. **Munging**:
   - Always munge summary statistics before analysis
   - Use `--merge-alleles` with HapMap3 SNP list

5. **Interpretation**:
   - Check regression intercept (should be close to 1)
   - Large intercepts suggest population stratification
   - Negative heritability estimates suggest poor quality data

## Common Issues and Solutions

### Issue: "WARNING: number of SNPs less than 200k"
**Solution**: Ensure summary statistics contain sufficient SNPs. Munge with HapMap3 SNP list.

### Issue: Negative heritability estimates
**Solution**: Check data quality, sample size reporting, and consider filtering low-quality SNPs.

### Issue: High LD Score regression intercept
**Solution**: May indicate population stratification or other confounding. Check genomic control lambda.

### Issue: "Could not find signed summary stats"
**Solution**: Ensure summary statistics include Z-scores, betas, or odds ratios.

## Performance Notes

- Heritability estimation: ~1-5 minutes for typical GWAS
- Genetic correlation: ~2-10 minutes per trait pair
- Partitioned heritability: ~10-30 minutes depending on annotations
- LD score computation: Varies with sample size and SNP count

## Building the Image

From this directory:

```bash
docker build -t biopsyk/ldsc:latest .
```

## References

### Primary Citations

If you use LDSC, please cite:

1. **LD Score Regression** (core method):
   - Bulik-Sullivan et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. *Nature Genetics*, 47, 291-295. [doi:10.1038/ng.3211](https://doi.org/10.1038/ng.3211)

2. **Genetic Correlation**:
   - Bulik-Sullivan et al. (2015). An atlas of genetic correlations across human diseases and traits. *Nature Genetics*, 47, 1236-1241. [doi:10.1038/ng.3406](https://doi.org/10.1038/ng.3406)

3. **Partitioned Heritability**:
   - Finucane et al. (2015). Partitioning heritability by functional annotation using genome-wide association summary statistics. *Nature Genetics*, 47, 1228-1235. [doi:10.1038/ng.3404](https://doi.org/10.1038/ng.3404)

4. **Stratified Heritability** (continuous annotations):
   - Gazal et al. (2017). Linkage disequilibrium-dependent architecture of human complex traits shows action of negative selection. *Nature Genetics*, 49, 1421-1427. [doi:10.1038/ng.3954](https://doi.org/10.1038/ng.3954)

### Additional Resources

- [Original LDSC Repository](https://github.com/bulik/ldsc)
- [Python 3 Version Repository](https://github.com/abdenlab/ldsc-python3)
- [LDSC Wiki](https://github.com/bulik/ldsc/wiki)
- [LD Hub](http://ldsc.broadinstitute.org/)

## Notes

- This image runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- Reference LD scores must be downloaded separately (not included in image)
- The image uses the Python 3 compatible version (v2.0.0) from abdenlab/ldsc-python3
- This is more modern than the original Python 2.7 version

## License

LDSC is licensed under GNU GPL v3.

## Authors

Original LDSC authors:
- Brendan Bulik-Sullivan (Broad Institute of MIT and Harvard)
- Hilary Finucane (MIT Department of Mathematics)

Python 3 port contributors:
- Thomas Reimonn and the abdenlab team
