# ReadStat Docker Image

This Docker image provides ReadStat, a command-line tool for reading and writing data files from statistical packages like SAS, Stata, and SPSS.

## Description

ReadStat is a C library for reading and writing data files from statistical packages. It allows you to convert between different statistical data formats (.dta, .por, .sav, .sas7bdat) and export to common formats like CSV and Excel.

## Version Information

- Internal image version: 1.0.0
- ReadStat: Latest stable version from GitHub (can be verified by running `readstat --help`)

## Usage

### Basic Usage

```bash
# Convert Stata file to CSV
docker run -v $(pwd):/data biopsyk/readstat readstat input.dta output.csv

# Convert SPSS file to Excel
docker run -v $(pwd):/data biopsyk/readstat readstat input.sav output.xlsx

# Convert SAS file to Stata format
docker run -v $(pwd):/data biopsyk/readstat readstat input.sas7bdat output.dta
```

### Using the convenience wrapper

```bash
# Use the readstat-docker wrapper for cleaner syntax
docker run -v $(pwd):/data biopsyk/readstat readstat-docker input.dta output.csv
```

### Using with Singularity

```bash
# Pull the image
singularity pull docker://biopsyk/readstat:1.0.0

# Convert files
singularity exec readstat_1.0.0.sif readstat input.dta output.csv
```

## Supported Formats

### Input Formats
- `.dta` - Stata files
- `.por` - SPSS portable files  
- `.sav` - SPSS system files
- `.sas7bdat` - SAS data files

### Output Formats
- `.dta` - Stata files
- `.por` - SPSS portable files
- `.sav` - SPSS system files  
- `.sas7bdat` - SAS data files
- `.xlsx` - Excel files
- `.csv` - Comma-separated values

## Examples

```bash
# Convert multiple files in a directory
for file in *.dta; do
    docker run -v $(pwd):/data biopsyk/readstat readstat "$file" "${file%.dta}.csv"
done

# Check file information without conversion
docker run -v $(pwd):/data biopsyk/readstat readstat --help
```

## Notes

- The container runs as a non-root user for security
- The working directory is set to `/data`
- Input files should be mounted to `/data` in the container
- The container includes only the necessary runtime dependencies to minimize size
- File paths should be relative to the mounted directory

## Building the Image

From this directory:

```bash
docker build -t biopsyk/readstat:latest .
```

## References

- [ReadStat GitHub Repository](https://github.com/WizardMac/ReadStat)
- [Original Docker implementation](https://github.com/jbn/readstat) 