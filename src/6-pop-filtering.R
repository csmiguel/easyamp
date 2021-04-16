###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: reformat and filter genotypes from population genetics
#PROJECT: amplicon-genotyping
###.............................................................................
library(dplyr)
library(pegas)

#load genotypes
gen_genind <- readRDS("data/intermediate/genotypes-genind.rds")

#load filtering parameters
source("src/parameters/pop-filtering.r")
#load filtering functions
source("src/functions/pop-filt.r")
#read ploidy
# I do a trick to prevent populating the working env with extra variables
get_ploidy <- function() {
  source("src/parameters/parameters.r", local = T)
  ploidy
}
ploidy <- get_ploidy()

sink("output/log.txt", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on src/5-pop-filtering.R\n")
#start filtering
pop_filtered_genotypes <-
  gen_genind %>%
  # filter out individuals according to their missing data
    filter_ind_missingness() %>%
  # filter out loci according to their missing data
    missingness_locus() %>%
  # remove monorphic loci
    remove_monomorphic() %>%
  # filter out loci not in HW equilibrium
    nonhw()
sink()

#export to structure
hierfstat::genind2hierfstat(pop_filtered_genotypes,
                            pop = adegenet::indNames(pop_filtered_genotypes)) %>%
  hierfstat::write.struct(fname = "output/output5-structure.str",
                          ilab = indNames(pop_filtered_genotypes),
                          MARKERNAMES = T,
                          pop = pop(pop_filtered_genotypes))
#save filtered genotypes
saveRDS(pop_filtered_genotypes, "data/output/output6-pop-genotypes-genind.rds")
