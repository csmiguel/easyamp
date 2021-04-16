[![DOI](https://zenodo.org/badge/279770729.svg)](https://zenodo.org/record/3945855)
# EasyAmp: a simple pipeline to genotype multilocus amplicon data from High Throughput Sequencing
The goal of this repository is to provide an easy tool to genotype loci from amplicon sequencing coming from High-Throughput Sequencing (HTS), such as from Illumina MiSeq.
* input: R1 and R1 amplicon reads from HTS demultiplexed by sample + primer sequences.
* output: individual x loci genotype matrix +  alleles in FASTA + read count per allele (and more).

It provides a modular and customizable workflow. The core of the genotyping module relies on DADA2, an R package designed for determining Amplicon Sequence Variants (ASVs) in metabarcoding.

![Alt text](etc/flow.svg?raw=true&sanitize=true)

## what it does:
* it genotypes amplicon reads from diploid data.
* it takes as input overlapping and non-overlapping reads.
* it can genotype simultaneously samples from the same or different species.
* various output formats: STRUCTURE, genotype table and FASTA.

## what it does **not** do:
* current version only works with **paired-end** reads.
* not for shotgun data.
* current version only works on **diploid** and in single-copy genes.
* not for genotyping gene families (e.g. MHC).
* most time consuming processes run in parallele, but others do not.

## set up

Clone the current repository.
You should have installed and added to the $PATH the following programs (also available in the `src` folder):
[cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) version 2.1 or above, and
`mmv` (MAC `brew install mmv`; LINUX `sudo apt-get install mmv`)

The scripts are dependent on the following R packages with their dependencies: `dada2 seqinr dplyr assertthat adegenet tibble reshape2 xlsx plyr hierfstat tidyr magrittr pegas`.
For more info see [R Session Info](`etc/R-session-info.txt`)

Raw sequences consist of demultiplexed reads (one file per sample). Move your R1 and R2 reads to `data/raw`.
Provide information on primers and loci as plain text in `data/raw/primers.csv`. Formatting instructions are written within each file.

All names of FASTQ files used should meet a given format. Please, follow instructions and edit the file `data/raw/id-match`. This file will be used to bulk rename all FASTQ files to meet the required formatting using `mmv`.

The names of loci and samples **must** only contain alphanumeric characters.

## repository structure

`src` contains scripts and source files.
`output`is an empty folder where the output of the scripts will be saved.
`data/intermediate` is where intermediate files from analysis are stored.
`data/raw` is where raw sequences and templates are stored.

## running the scripts

Scripts should be run in order, starting from 0 (`src/0-checks.R`) in the Terminal/Console, or inside R if applies, and using relative paths.

`src/0-checks.R` confirms you have all the installed dependencies and creates input files with primers for cutadapt.

`src/1-trim-reads.sh` renames FASTQ to meet input format and runs cutadapt. It demultiplexes loci from R1 and R2 files into separate FASTQ files in `data/intermediate`. Log from cutadapt can be found in `output/cutadapt.txt`.

`src/2-genotyping.R` does the genotyping for each locus across all samples. Reads with ambiguities are removed. If reads do not overlap (i.e. the amplicon is longer than  ~R1 + R2), then they are merged by adding 10 N in between R1 and R2 reads. Final genotype calls are given using thresholds on minimium coverage and allele balance `src/parameters/parameters.r`.
heterozygous (AB): two alleles found passing thresholds in `src/parameters/parameters.r`.
homozygous (AA): one allele found. The read count for that allele minus the `cov` threshold in `src/parameters/parameters.r` is above 5.
hemizygous (A-): one allele found. The read count for that allele minus the `cov` threshold in `src/parameters/parameters.r` is below 5.

It creates:
* `output/output1-genotypes.tsv`, with genotypes as a table.
* `output/output2-allele-freqs.txt`, table with samples, loci and read count for each allele.
* `output/output3-alleles.fasta`, fasta with all different alleles.
* `output/output4-individuals-x-alleles.fasta`, fasta with all alleles per individual.

YOU CAN STOP HERE. However, if you have populational data you might want to run scripts below for further filtering of the genotypes.

`src/3-filtering.R` filters the genotypes. It drops individuals or loci with missing data above/below thresholds in `src/parameters/parameters.r`. It removes monomorphic loci and loci not in HWE. Edit filtering thresholds in `src/parameters/parameters.r`.
A report of the filtering is created in `output/population-filtering.log`.

`src/4-reformat.R` reformats population genetics data. It requires having run `src/3-filtering.R`. It generates:
* `output/output5-filtered-genotypes.tsv`, a table with the genotypes.
* `output/output6-genotypes.str`, input genotypes for STRUCTURE.
* `output/output7-genind.rds`, genotypes as a `genind` object.


## some notes

Any time, you can delete intermediate results `rm -rf data/intermediate` and output `rm -rf output`, to re-do analysis from raw data.
