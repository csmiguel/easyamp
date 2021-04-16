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
#loci
loci <- readLines("data/intermediate/loci")
#vector to name of alleles
al_names_all <- as.character(101:999)
#log genotyping
zz <- file("output/log.txt", open = "a")
sink(zz, type = "output", append = T)
sink(zz, type = "message", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on 3-genotyping.R\n")

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

#bind genotype data
gen <-
  genotypes_notnull %>%
  plyr::ldply(function(x) x$asvs, .id = NULL) %>% tibble::as_tibble()
#bind allele data
al <-
  genotypes_notnull %>%
  plyr::ldply(function(x) x$al_seq, .id = NULL) %>%
  dplyr::rename(allele = name) %>% tibble::as_tibble()
# create list
genotypes_notnull <- list(asvs = gen, al_seq = al)

#save objects
saveRDS(genotypes_notnull, file = "data/intermediate/genotypes_notnull.rds")
