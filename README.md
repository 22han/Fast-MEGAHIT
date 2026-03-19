# Fast-MEGAHIT

Fast metagenomic assembler with kmerlight integration for efficient k-mer counting.

## Quick Start

### 1. Clone and Build
```bash
git clone https://github.com/22han/Fast-MEGAHIT.git
cd Fast-MEGAHIT
```
### 2. Run the Pipeline with One Command
```bash
./parallel.sh /path/to/reads_1.fastq /path/to/reads_2.fastq /path/to/output_dir
```
reads_1.fastq, reads_2.fastq – forward and reverse reads (FASTQ format)

output_dir – directory where all results will be saved


### 3.Requirements
CMake ≥ 3.10

Python ≥ 3.9

GCC/Clang compiler (with C++11 support)

Memory: At least 8 GB for small test datasets. For large metagenomic assemblies (e.g., tens of millions of reads), hundreds of GB of RAM may be required.

CPU: Multi-core recommended; the pipeline scales well with many threads (tested with up to 144 threads).

Example Dataset (Optional)
To test the pipeline with a small public dataset, you can download ERR14952688 from ENA:

bash
### Install fastq-dump if needed (from sra-tools)
```
fastq-dump --split-files ERR14952688
```
Then run:

```bash
./parallel.sh ERR14952688_1.fastq ERR14952688_2.fastq ./test_output
```
