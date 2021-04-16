###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: applies dada2::FilterAndTrim and outputs reads in intermediate/filtered
#PROJECT: amplicon-genotyping
###.............................................................................
library(dada2)
library(dplyr)

#load parameters for genotype calling
source("src/parameters/parameters.r")
#load function to cleas spurious ASVs
source("src/functions/trimandfilter.r")
#loci
loci <- readLines("data/intermediate/loci")

#"data" path for filtered sequences downstream
filt_path <- "data/intermediate/filtered"
if (!file_test("-d", filt_path)) dir.create(filt_path)

#log genotyping
zz <- file("output/log.txt", open = "a")
sink(zz, type = "output", append = T)
sink(zz, type = "message", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on src/2-genotyping.R\n")

# filter and trim reads
filtering_output <- loci %>%
  lapply(function(locus) {
    if (sepe == "pe") h <- trim_pe(locus)
    if (sepe == "se") h <- trim_se(locus)
    return(h)
  }) %>% setNames(loci)

cat("\nFiltered reads after truncation (dada2::TrimAndFilter)\n")
filtering_output

closeAllConnections()
