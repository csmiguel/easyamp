###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: reformat and filter genotypes from population genetics data
#PROJECT: amplicon-genotyping
###.............................................................................
library(dplyr)

#load genotypes
genotypes <- readRDS("data/intermediate/genotypes.rds")
init_gen <- genotypes
#load filtering parameters
source("src/parameters/parameters.r")
#load filtering functions
source("src/functions/gen-filt.r")

#remove non-scored loci
unscored_loci <- unscored(genotypes)
  #   apply filter
  genotypes <- genotypes[, !unscored_loci]

#filter out list of individuals
  if (is.character(ind_to_remove)){
    genotypes <- genotypes[!(rownames(genotypes) %in% ind_to_remove), ]
  }
  
#filter out individuals according to their missing data
ind_miss <- ind(genotypes)
  # apply filter
  genotypes <-
    genotypes[ind_miss, ]

#filter out loci according to their missing data
locus_miss <- locus(genotypes)
  # apply filter
  genotypes <-
    genotypes[, locus_miss]

  
#remove monorphic loci (applied the latest)
monomorphic_loci <- monomorphic(genotypes)
# apply filter
genotypes <- genotypes[, !monomorphic_loci]

#filter out loci non in HW equilibrium
hw_rm <- nonhw(genotypes, pvalue_hw)
# apply filter
genotypes <-
  genotypes[, names(genotypes) %in% hw_rm]

#create report
sink("output/population-filtering.log")
  cat("input file has", nrow(init_gen), "individuals and ",
      ncol(init_gen), "loci\n")
  cat(ifelse(all(!unscored_loci), "0",
             paste(names(unscored_loci)[unscored_loci], collapse = " ")),
      "has/have not been genotyped and has been removed", "\n")
  cat(ifelse(all(ind_miss), "0",
           paste(names(ind_miss)[!ind_miss], collapse = " ")),
    "samples has/have missing data above threshold",
    max_missing_ind, "\n")
  cat(ifelse(all(locus_miss), "0",
           paste(names(locus_miss)[!locus_miss], collapse = " ")),
    "loci has/have missing data above threshold",
    max_missing_loci, "\n")
  cat(ifelse(all(!monomorphic_loci), "0",
             paste(names(monomorphic_loci)[monomorphic_loci], collapse = " ")),
      "loci is/are monomorphic and have been removed", "\n")
  cat(paste(hw_rm, collapse = " "),
      "has/have been remved because it/they was/were not in HW equilibrium at p =", pvalue_hw,"\n")
  cat("output file has", nrow(genotypes), "individuals and ",
      ncol(genotypes), "loci", "\n")
sink()

#save filtered genotypes
saveRDS(genotypes, "data/intermediate/genotypes-filtered.rds")
