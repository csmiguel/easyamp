###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: reformat genotypes
#PROJECT: amplicon-genotyping
###.............................................................................
library(adegenet)
library(hierfstat)
library(dplyr)
library(magrittr)
#load reformatting options
source("src/parameters/parameters.r")
#read genotypes
genotypes <-
  readRDS("data/intermediate/genotypes-filtered.rds") %>%
  apply(1, function(x) { #replace calls with one allele only, [a-z]0, with [a-z]
  position <- grep("[a-z]0", x)
  x[position] %<>% stringr::str_remove("0")
  x
}) %>% t()

#1. convert to genind object
gen_ad <- adegenet::df2genind(
            genotypes,
            ncode = 2,
            ind.names = row.names(genotypes),
            loc.names = names(genotypes),
            sep = "",
            NA.char = "0",
            ploidy = 2,
            type = "codom")
saveRDS(gen_ad, "output/output7-genind.rds")

#2. export in STRUCTURE format
# create a second gen_ad object
gen_ad2 <- gen_ad
  # hierfstat requires names of alleles to be be integers:
gen_ad2@all.names %<>%
  sapply(function(x) {
    plyr::mapvalues(
      x,
      letters,
      1:length(letters))
  })
# export in STRUCTURE format
hierfstat::genind2hierfstat(gen_ad2,
                            pop = adegenet::indNames(gen_ad2)) %>%
  hierfstat::write.struct(fname = "output/temp")
# replace loci with only one allele genotyped with missing data
system("cat output/temp | sed 's/ 0 / -9 /g' > output/output6-genotypes.str; rm output/temp")

#3. write to table
write.table(genotypes, file = "output/output5-filtered-genotypes.tsv", quote = F)
