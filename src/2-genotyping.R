###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: dada2
#PROJECT: amplicon-genotyping
###.............................................................................
library(dada2)
library(dplyr)
library(ape)
library(kmer)
library(magrittr)

#load parameters for genotype calling
source("src/parameters/parameters.r")
#load function to cleas spurious ASVs
source("src/functions/genotyping.r")

#samples
all_samples <- readLines("data/intermediate/samples-list")
#loci
loci <- readLines("data/intermediate/loci")
#vector to name of alleles
al_names_all <- as.character(101:999)
#create path for filtered sequences downstream
filt_path <- "data/intermediate/filtered"
if (!file_test("-d", filt_path)) dir.create(filt_path)

#log genotyping
zz <- file("output/log.txt", open = "a")
sink(zz, type = "output", append = T)
sink(zz, type = "message", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on src/2-genotyping.R\n")

#genotyping. See function details for further info.
genotypes <-
  loci %>%
  lapply(function(locus) genotype(locus)) %>%
    setNames(loci)
#remove NULL loci; ie those that have not been genotyped
notnull <- genotypes %>% sapply(function(x) !is.null(x))
notnull_loci <- loci[notnull]
genotypes_notnull <- genotypes[notnull]
#print loci removed
cat("\nLoci removed:")
print(loci[!(loci %in% notnull_loci)])
sink()

#save objects
saveRDS(genotypes_notnull, file = "data/intermediate/genotypes_notnull.rds")
