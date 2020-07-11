###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: run checks on files necessary and create files
#PROJECT: amplicon-genotyping
###.............................................................................
library(asserthat)
library(dplyr)

#are all necessary packages installed?
h <- installed.packages()
assertthat::assert_that(
  !(match(
    c("dada2", "seqinr", "dplyr", "assertthat", "adegenet", "tibble", "reshape2", "xlsx", "plyr", "hierfstat"), h[, 1], nomatch = NA) %>%
    sum() %>% is.na()),
  msg = "you need to install all necessary R packages"
)
#1. PRIMERS: checks on input data and create fasta for forward and reverse primers
#read data
path_txt <- "data/raw/data"
path_excel <- "data/raw/data.xls"
if (file.exists(path_txt)) {
  data <- read.csv(path_txt)
  } else if (file.exists(path_excel)) {
  data <- xlsx::read.xlsx(path_excel, 1)
  } else {
  warning("data/raw/data or data/raw/data.xls not found. Provide primer sequences and loci names in those files")
  }

#assert data format looks ok
assertthat::assert_that(
  !is.null(data) & apply(data, 2, length) %>%
                    unique() %>%
                    length() == 1,
                    msg = paste("check", path_txt, "or", path_excel))

#create fasta file with forward primers
data %>%
  dplyr::select(1) %>%
  dplyr::pull() %>%
  lapply(function(x) x) %>%
  seqinr::write.fasta(names = data[, 3], file.out = "data/intermediate/forward")
#create fasta file with forward primers
data %>%
  dplyr::select(2) %>%
  dplyr::pull() %>%
  lapply(function(x) x) %>%
  seqinr::write.fasta(names = data[, 3], file.out = "data/intermediate/reverse")
#create text file with names of loci
write(data$locus, file = "data/intermediate/loci")

#2. raw sequences
#raw sequences should be placed in data/raw/ and be in fastq format
#compression is optional
#fill in the file raw/data/id-match
# are there any fastq files present in data/raw?

assertthat::assert_that(dir("data/raw", "fastq") %>% length() > 1,
        msg = "no fastq files found on data/raw")

#has raw/data/id-match be filled?

assertthat::assert_that(
  readLines("data/raw/id-match") %>%
      grep(pattern = "^data/raw") %>%
      length() > 0,
  msg =  "nothing to rename found on data/raw/id-match")
