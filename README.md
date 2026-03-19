

# Fast-MEGAHIT

Fast metagenomic assembler with kmerlight integration for efficient k-mer counting.

 ## Quick Start

### 1. Clone Repository
```
git clone https://github.com/22han/Fast-MEGAHIT.git 
cd Fast-MEGAHIT
```
### 2. Build from Source
```
make
```
### 3. Quality Control (Preprocessing)
```
python3 /path/to/pre.py SRR33980743_1.fastq SRR33980743_2.fastq
```
### 4. Parallel Processing with kmerlight
```
./parallel.sh \
    kmerlight-master \
    parallel_analysis_output \
    clean_1.fastq \
    clean_2.fastq
```
### 5. Run Assembly
```
 ./megahit \
    -1 clean_1.fastq \
    -2 clean_2_rc.fastq \
    -o assembly_output \
    --k-list 29,31,33,51,53,77,87,103,117,135 \
    --min-contig-len 2000
```
## MEGAHIT Options
-1, -2: Forward/reverse read files (FASTQ format)

-o: Output directory

--k-list: K-mer sizes (comma-separated, odd numbers)

--min-contig-len: Minimum contig length (default: 2000)

## Example Dataset
Use sample.fastq for testing and validation.

## Requirements
CMake >= 3.10
Python >= 3.6
GCC/Clang compiler
8GB+ RAM (for large datasets)
