# Multilocus genotyping from amplicon sequencing
The goal of this repository is to provide an easy tool to genotype loci from amplicon sequencing coming from High-Throughput Sequencing (HTS), such as from Illumina MiSeq.
* Input: amplicon reads from HTS demultiplexed by sample. Loci primer sequences.
* Output: individual x loci genotype matrix, plus alleles in FASTA.

It provides a modular and customizable workflow. The core of the genotyping module relies on dada2, an R package designed for determining Amplicon Sequence Variants (ASVs) in metabarcoding.

## what it does:
* it genotypes amplicon reads from diploid data.
* it takes as input overlapping and non-overlapping reads.
* it can genotype simultaneously samples from the same or different species.
* various output format: STRUCTURE, genotype table and FASTA.

## what it does **not** do:
* current version only works with **paired-end** reads.
* not for shotgun data.
* current version only works on **diploid** and in single-copy genes.
* do not use for genotyping gene families (eg MHC)
* most time consuming processes run in parallele, but others do not.

## set up

Clone the current repository.
You should have installed and added to the $PATH the following programs (also available in the `src` folder):
[cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) version 2.1 or above
mmv (MAC `brew install mmv`; LINUX `sudo apt-get install mmv')

The scripts are dependent on the following R packages with their dependencies: `dada2 seqinr dplyr assertthat adegenet tibble reshape2 xlsx plyr hierfstat tidyr magrittr`.
For more info see [R Session Info](`etc/sessionInfo.txt`)

Raw sequences consist of demultiplexed reads (one file per sample). Move your R1 and R2 reads to `data/raw`.
Provide information on primers and loci in the Excel file `data/raw/data.xls`, or in the as plain text in `data/raw/data`. Formatting instructions are written within each file.

All names of fastq files used should meet a given format. Please, follow instructions and edit the file `data/raw/id-match`. This file will be used to bulk rename all fastq files to meet the required formatting using `mmv`.

## repository structure

`src` contains scripts and source files.
`output`is an empty folder where the output of the scripts will be saved.
`data/intermediate` is an empty folder where intermediate files are stored.
`data/raw` is where raw sequences and templates are stored.

## running the scripts

Scripts should be run in order from 0 (`src/0-checks.R`)to 5 (`code/05-fasfda.R`) in the Terminal/Console, or inside R if applies, using relative paths.

`src/0-checks.R` confirms you have all the installed dependencies and creates input files with primers for cutadapt.

`src/1-trim-reads.sh` remanes fastq to meet input format and runs cutadapt. It demultiplexes loci from R1 and R2 files into separate fastq files in `data/intermediate`. Log from cutadapt can be found in `output/cutadapt.txt`.

`src/2-genotyping.R` does the genotyping for each locus across all samples. Reads with ambiguities are removed. If reads do not overlap (i.e. the amplicon is longer than  ~R1 + R2), then they are merged by adding 10 N in between R1 and R2 reads. Final genotype calls are given using thresholds on minimium coverage and allele balance `src/parameters/parameters.r`.
heterozygous (AB): two alleles found passing tresholds in `src/parameters/parameters.r`.
homozygous (AA): one allele found. The read count for that allele minus the `cov` threshold in `src/parameters/parameters.r` is above 5.
hemizygous (A-): one allele found. The read count for that allele minus the `cov` threshold in `src/parameters/parameters.r` is below 5.

It creates:
* `output/alleles.fasta` fasta with all different alleles.
* `output/allele-freqs.txt` table with samples, loci and read count for each allele.
* `output/individuals-x-alleles.fasta` fasta with all alleles per individual.
* `output/genotypes.tsv` with genotypes as a table.

YOU CAN STOP HERE. However, if you have populational data you might want to run scripts below to filter genotypes given their missing data and write output for STRUCTURE 

`src/3-filtering.R` applies filters to the genotypes. It drops individuals or loci with missing data above/below thresholds in `src/parameters/parameters.r`. It removes monomorphic loci. Edit filtering thresholds in `src/parameters/parameters.r`.
A report of the filtering is created in `output/population-filtering.log`.

`src/4-reformat.R` reformats population genetics data. It requires running `src/3-filtering.R`. It generates:
* `output/filtered-genotypes.tsv`, a table with the genotypes.
* `output/structure.str`, input genotypes for STRUCTURE.
* `output/genind.rds`, genotypes as a `genind` object.


## some notes

At any time,  you can delete intermediate results `rm -rf data/intermediate` and output `rm -rf output`, to re-do analysis from raw data.
