# SOLAR-Eclipse Docker Image

This Docker image provides SOLAR-Eclipse, a comprehensive software package for genetic variance components analysis and pedigree-based genetic studies.

## Description

SOLAR-Eclipse (Sequential Oligogenic Linkage Analysis Routines) is an extensive, flexible software package for genetic variance components analysis. It extends the capabilities of the original SOLAR software to support modern neuroimaging and genetic studies, including mega and meta-analyses of genetic imaging data.

## Version Information

- Internal image version: 1.0.0
- SOLAR-Eclipse version: 9.0.0 (Static Linux)
- Built from: https://www.nitrc.org/projects/se_linux/

## Key Features

- **Variance Components Analysis**: Decompose phenotypic variance into genetic and environmental components
- **Heritability Estimation**: Calculate narrow-sense and broad-sense heritability
- **Linkage Analysis**: Multipoint and marker-specific IBD calculations for pedigrees of any size
- **Association Analysis**: SNP-based QTN (Quantitative Trait Nucleotide) and QTLD (Quantitative Trait Locus Disequilibrium) tests
- **Polygenic Modeling**: Account for family structure and random genetic effects
- **Covariate Screening**: Systematic evaluation of potential confounders
- **Imaging Genetics**: Optimized for neuroimaging phenotypes (voxel-wise, ROI-based)
- **Complex Pedigrees**: Handles pedigrees of arbitrary size and complexity
- **Multi-trait Analysis**: Analyze multiple phenotypes simultaneously
- **Oligogenic Analysis**: Model multiple genetic loci
- **Environmental Effects**: Include household, dominance, and interaction effects

## Usage

### Basic Commands

```bash
# Show help
docker run biopsyk/solar

# Run SOLAR interactively
docker run -it -v $(pwd):/data biopsyk/solar solar

# Run SOLAR script
docker run -v $(pwd):/data biopsyk/solar solar < analysis_script.tcl
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/solar:1.0.0

# Run SOLAR interactively
singularity exec solar_1.0.0.sif solar

# Run SOLAR script
singularity exec solar_1.0.0.sif solar < analysis_script.tcl
```

## Common Workflows

### 1. Basic Heritability Estimation

Calculate the heritability of a quantitative trait in a pedigree:

```tcl
# Load pedigree data
load pedigree myped.csv

# Load phenotype data
load phenotypes phenotypes.csv

# Select the trait to analyze
trait height

# Run polygenic analysis
polygenic

# Display results
phen

# Save results
outfile height_h2.out
phen
outfile /dev/tty
```

Docker command:
```bash
docker run -v $(pwd):/data biopsyk/solar solar < heritability.tcl
```

### 2. Covariate Screening

Identify significant covariates for a trait:

```tcl
load pedigree myped.csv
load phenotypes phenotypes.csv

trait height

# Screen all covariates
covariate age sex bmi
polygenic -screen

# Use selected covariates
covariate age sex
polygenic

# Display final model
model
```

### 3. Bivariate Genetic Correlation

Estimate genetic correlation between two traits:

```tcl
load pedigree myped.csv
load phenotypes phenotypes.csv

# Analyze trait 1
trait height
polygenic
phen

# Analyze trait 2
trait weight
polygenic
phen

# Genetic correlation analysis
rhog height weight
```

### 4. Linkage Analysis

Perform genome-wide linkage scan:

```tcl
load pedigree myped.csv
load phenotypes phenotypes.csv
load markers markers.csv
load map genetic_map.csv

trait height
covariate age sex

# Calculate multipoint IBD matrices
multipoint -b -i

# Run linkage analysis
polygenic
multipoint

# Extract LOD scores
lodsummary lod_scores.txt
```

### 5. SNP Association Analysis (QTLD)

Test for association with SNP genotypes:

```tcl
load pedigree myped.csv
load phenotypes phenotypes.csv
load snp_data snps.csv

trait height
covariate age sex

# Run polygenic model
polygenic

# Test each SNP
for {set i 1} {$i <= $nsnps} {incr i} {
    snp $i
    assoc
}

# Or test a specific SNP
snp rs1234567
assoc
```

### 6. Imaging Genetics Analysis

Analyze neuroimaging phenotypes with genetic data:

```tcl
load pedigree family.csv
load phenotypes imaging_phenotypes.csv

# Voxel-wise analysis example
define voxel_list = [list vox_1 vox_2 vox_3 ... vox_n]

foreach voxel $voxel_list {
    trait $voxel
    covariate age sex scanner
    polygenic
    phen
}

# ROI-based analysis
trait cortical_thickness_L_frontal
covariate age sex intracranial_volume
polygenic
phen
```

### 7. Mega-Analysis of Multiple Sites

Combine data from multiple cohorts:

```tcl
# Load data from multiple sites
load pedigree site1.ped
load phenotypes site1.phe

load pedigree site2.ped
load phenotypes site2.phe

# Account for site effects
trait phenotype
covariate age sex site
polygenic

# Test for site heterogeneity
house site
polygenic
```

## SOLAR Command Reference

### Data Loading Commands

```tcl
load pedigree <file>     # Load pedigree structure (ID, FA, MO, SEX)
load phenotypes <file>   # Load phenotype data
load markers <file>      # Load marker genotypes
load map <file>          # Load genetic map
load snp_data <file>     # Load SNP data
```

### Analysis Setup

```tcl
trait <name>             # Select trait for analysis
covariate <vars>         # Specify covariates
house <var>              # Specify household effect
```

### Core Analysis Commands

```tcl
polygenic                # Run polygenic model
polygenic -screen        # Screen covariates
assoc                    # SNP association test
multipoint               # Linkage analysis
multipoint -i            # Calculate IBD matrices
rhog <trait1> <trait2>   # Genetic correlation
```

### Model Specification

```tcl
model                    # Display current model
define <var> = <value>   # Define SOLAR variable
option <name> <value>    # Set SOLAR option
```

### Output Commands

```tcl
phen                     # Display phenotypic model results
outfile <file>           # Redirect output to file
outfile /dev/tty         # Redirect output to terminal
lodsummary <file>        # Write LOD score summary
```

### Utility Commands

```tcl
help                     # Display help
help <command>           # Help for specific command
pwd                      # Print working directory
cd <dir>                 # Change directory
ls                       # List directory contents
quit                     # Exit SOLAR
```

## Input File Formats

### Pedigree File Format

Standard pedigree format with tab or comma-separated values:

```
ID      FA      MO      SEX
1       0       0       1
2       0       0       2
3       1       2       1
4       1       2       2
```

Required columns:
- **ID**: Individual identifier
- **FA**: Father ID (0 if founder)
- **MO**: Mother ID (0 if founder)
- **SEX**: Sex (1=male, 2=female)

### Phenotype File Format

Tab or comma-separated file with phenotypes:

```
ID      height  weight  age     sex
1       175.5   70.2    45      1
2       162.3   58.4    43      2
3       180.1   75.8    20      1
4       168.7   62.1    18      2
```

First column must be ID matching the pedigree file.

### Marker/Map Files

Marker genotype file:
```
ID      SNP1    SNP2    SNP3
1       1/1     1/2     2/2
2       1/2     2/2     1/2
```

Genetic map file:
```
chromosome  position_cM  marker_name
1           0.0          SNP1
1           1.5          SNP2
1           3.2          SNP3
```

## Output Interpretation

### Heritability Results

Key output values from `phen` command:

- **H2r**: Narrow-sense heritability (additive genetic variance / total variance)
- **H2r SE**: Standard error of heritability estimate
- **P-value**: Significance test (H2r > 0)
- **Residual SD**: Standard deviation of environmental variance
- **LogLike**: Log-likelihood of the model

Example output:
```
H2r = 0.68 ± 0.12  (p = 2.3e-08)
```
Interpretation: 68% of trait variance is due to additive genetic factors.

### Linkage Analysis Results

- **LOD score**: Log-odds score for linkage
  - LOD > 3.3: Significant linkage (genome-wide)
  - LOD > 2.0: Suggestive linkage
  - LOD > 1.0: Nominal linkage
- **Position**: Genetic map position (cM)
- **Chromosome**: Chromosome number

### Association Analysis Results

- **Chi-square**: Test statistic for association
- **P-value**: Significance of association
- **Effect size**: Estimated effect of SNP on trait

## Best Practices

### 1. Data Preparation

- **Quality Control**: Remove individuals/markers with >10% missing data
- **Outliers**: Identify and investigate extreme phenotype values
- **Pedigree Errors**: Verify pedigree structure with genetic data
- **Trait Distribution**: Transform non-normal traits (log, inverse normal)

### 2. Covariate Selection

- **Always Include**: Age, sex as appropriate
- **Screening**: Use `-screen` option to identify significant covariates
- **Multicollinearity**: Avoid highly correlated covariates
- **Stratification**: Include population principal components if needed

### 3. Model Building

- **Start Simple**: Begin with basic polygenic model
- **Add Complexity**: Progressively add covariates, dominance, etc.
- **Compare Models**: Use likelihood ratio tests
- **Check Convergence**: Ensure models converge properly

### 4. Imaging Genetics

- **Multiple Testing**: Correct for family-wise error rate or FDR
- **Covariates**: Include scanner, site, ICV for brain measures
- **Batch Effects**: Account for acquisition differences
- **Voxel Selection**: Pre-filter voxels by variance or heritability

### 5. Performance Optimization

- **Large Pedigrees**: SOLAR handles arbitrarily large pedigrees efficiently
- **Parallel Processing**: Run independent analyses in parallel
- **IBD Caching**: Reuse calculated IBD matrices when possible
- **Memory**: Static build included in this image minimizes memory usage

## Common Issues and Solutions

### Issue: "Trait has too many missing values"
**Solution**: Remove individuals with missing data or use SOLAR's missing data handling:
```tcl
option MissingPhenotype 999
```

### Issue: "Matrix not positive definite"
**Solution**: Pedigree structure may be problematic. Check for:
- Inbreeding loops
- Very large sibships
- Incorrect pedigree links

### Issue: Convergence problems
**Solution**:
- Check trait distribution (transform if needed)
- Remove extreme outliers
- Simplify model (reduce covariates)
- Verify pedigree structure

### Issue: Very low/high heritability estimates
**Solution**:
- Verify phenotype quality and units
- Check covariate specifications
- Ensure sufficient family structure in data
- Consider shared environmental effects (household)

## Performance Notes

- **Pedigree Size**: Handles pedigrees from families to large cohorts (10,000+ individuals)
- **Speed**: Polygenic analysis typically completes in seconds to minutes
- **Memory**: Static build minimizes memory requirements
- **Parallelization**: Multiple independent analyses can run simultaneously

## Example Complete Analysis

```tcl
#!/usr/bin/env solar

# Complete heritability analysis script

# Load data
load pedigree family_structure.csv
load phenotypes measurements.csv

# Setup
trait cortical_thickness
covariate age sex intracranial_volume

# Screen for additional covariates
covariate age sex intracranial_volume education bmi
polygenic -screen

# Final model with selected covariates
covariate age sex intracranial_volume
polygenic

# Output results
outfile results/cortical_thickness_h2.txt
phen
model
outfile /dev/tty

# Bivariate analysis with another trait
trait hippocampal_volume
covariate age sex intracranial_volume
polygenic
phen

# Genetic correlation
rhog cortical_thickness hippocampal_volume

quit
```

Run with Docker:
```bash
docker run -v $(pwd):/data biopsyk/solar solar < analysis_script.tcl
```

## Building the Image

From this directory:

```bash
docker build -t biopsyk/solar:latest .
```

## References

### Primary Citations

- Almasy, L., & Blangero, J. (1998). Multipoint quantitative-trait linkage analysis in general pedigrees. *American Journal of Human Genetics*, 62(5):1198-1211. [doi:10.1086/301844](https://doi.org/10.1086/301844)

- Kochunov, P., et al. (2019). SOLAR-Eclipse: A new neuroimaging data analysis package. *Frontiers in Neuroinformatics*, 13:16. [doi:10.3389/fninf.2019.00016](https://doi.org/10.3389/fninf.2019.00016)

### Additional Resources

- [SOLAR-Eclipse Official Website](https://solar-eclipse-genetics.org/)
- [NITRC Project Page](https://www.nitrc.org/projects/se_linux/)
- [NIH HPC SOLAR Manual](https://hpc.nih.gov/docs/solar-8.1.1/)
- [Solarius R Package](http://ugcd.github.io/solarius/) (R interface to SOLAR)

## Notes

- This image uses the static Linux build for maximum compatibility
- The image runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- SOLAR uses tcsh syntax for scripting
- Interactive mode is available with `-it` flags

## License

SOLAR-Eclipse is free, open-source software available for download and use in research. See the official website for licensing details.

## Authors

Original SOLAR developed by:
- Laura Almasy (University of Pennsylvania)
- John Blangero (Texas Biomedical Research Institute)

SOLAR-Eclipse extensions by:
- Peter Kochunov (University of Maryland)
- And contributors

## Support

For questions and support:
- Visit the [SOLAR-Eclipse website](https://solar-eclipse-genetics.org/)
- Check the [NIH HPC documentation](https://hpc.nih.gov/docs/solar-8.1.1/)
- Post issues on the NITRC project page
