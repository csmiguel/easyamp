###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: run checks on files necessary and create files
#PROJECT: amplicon-genotyping
###.............................................................................
library(assertthat)
library(dplyr)

# 1. install all necessary packages
h <- installed.packages()
pkg_needed <-
  c("dada2", "seqinr", "dplyr", "assertthat", "adegenet", "tibble",
    "reshape2", "xlsx", "plyr", "hierfstat", "kmer", "ape", "magrittr",
  "ggplot2", "gridExtra")
new_packages <- pkg_needed[!(pkg_needed %in% h[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
h <- installed.packages()

# 2. check there is a file with primers
#read data
path_csv <- "data/raw/primers.csv"

if (file.exists(path_csv)) {
  primers <- read.csv(path_csv, sep = ";")
  } else {
  warning("Provide data/raw/primers.csv sequences and loci names in those files")
  }

assertthat::assert_that(nrow(primers) > 0,
  msg = "fill data/raw/primers.csv with primer sequences and name of loci")
assertthat::assert_that(
  !is.null(primers) & apply(primers, 2, length) %>%
                    unique() %>%
                    length() == 1,
                    msg = paste("check", path_csv))
# 3. create fasta files with primers
# forward
primers %>%
  dplyr::select(1) %>%
  dplyr::pull() %>%
  lapply(function(x) x) %>%
  seqinr::write.fasta(names = primers[, 3], file.out = "data/intermediate/forward")
# reverse
primers %>%
  dplyr::select(2) %>%
  dplyr::pull() %>%
  lapply(function(x) x) %>%
  seqinr::write.fasta(names = primers[, 3], file.out = "data/intermediate/reverse")

#create text file with names of loci
write(as.character(primers$locus), file = "data/intermediate/loci")

#4. raw sequences
#raw sequences should be placed in data/raw/ and be in fastq format
#compression is optional
#fill in the file raw/data/id-match
# have fastq files been placed in data/raw?

assertthat::assert_that(dir("data/raw", "fastq") %>% length() > 1,
        msg = "no fastq files found on data/raw")

#has raw/data/id-match be filled?
assertthat::assert_that(
  readLines("data/raw/id-match.txt") %>%
      grep(pattern = "^data/raw") %>%
      length() > 0,
    msg =  "nothing to rename found on data/raw/id-match")

#log
hh <- pkg_needed[!(pkg_needed %in% h)]
logs <- "output/log.txt"
if (!file.exists(logs)) file.create(logs)
sink(logs)
cat(
  as.character(Sys.time()),
  "\nLOG on src/0-checks.R\n",
  paste(if (length(hh) == 0) "No" else hh,
   "package(s) need to be installed manually\n"),
   "data/raw/primers.csv exists and seems to have data formatted correctly\n",
   "fastq files seem to be present in data/raw\n",
   "fastq files seem to have a file data/raw/id-match.txt for renaming\n")
sink()
