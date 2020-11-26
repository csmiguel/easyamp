###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# November 2020
###.............................................................................
#GOAL: reformat genotypes and generate text outputs
#PROJECT: amplicon-genotyping
###.............................................................................
library(reshape2)
library(adegenet)
library(dplyr)

#read genotypes
genotypes <- readRDS("data/intermediate/genotypes_notnull.rds")

#read populations
popfile <- "data/raw/pops"
if (file.exists(popfile)) {
  populations <- read.csv(popfile,
    header = F, col.names = c("sample", "population"))
} else {
  populations <- NULL
  }
#read ploidy
# I do a trick to prevent populating the working env with extra variables
get_ploidy <- function() {
  source("src/parameters/parameters.r", local = T)
  ploidy
}
ploidy <- get_ploidy()
#bind genotype data
gen <-
  genotypes %>%
  plyr::ldply(function(x) x$asvs, .id = NULL) %>% tibble::as_tibble()
#bind allele data
al <-
  genotypes %>%
  plyr::ldply(function(x) x$al_seq, .id = NULL) %>%
  dplyr::rename(allele = name) %>% tibble::as_tibble()

#output 1: table with genotypes
#  |         | locus1 | locus2 |
#  |---------|--------|--------|
#  | sampleA | aa     | ab     |
#  | sampleB | bb     | ab     |

#reshape genotype data to:
h <- reshape2::acast(data = gen,
                formula = sample ~ locus + allele)
colnames(h) <- stringr::str_replace(colnames(h), "_", ".")

#create population vector
# assign a ficticious unique population to the data if no pop data is provided
if(is.null(populations)) {
  pops <- rep("1", nrow(h))
} else {
# if population data is provided, use it here:
  pops <-
    data.frame(sample = rownames(h)) %>%
    dplyr::left_join(populations, by = "sample") %>%
    dplyr::pull(population)
  }

#store genotypes in genind object
hh <- new("genind",
          tab = h,
          type = "codom",
          ploidy = ploidy,
          pop = pops)
#create hierfstat like dataframe
gen_df <-
  hierfstat::genind2hierfstat(hh) %>%
  dplyr::filter(pop != "dumpop")
#write output1
write.csv(gen_df, "output/output1-table-genotypes.csv")

#output2: table with samples, locus, allele and read count
# | sample  | locus  | allele | read count |
# |---------|--------|--------|------------|
# | sampleA | locus1 | a      | 1013       |
# | sampleA | locus1 | b      | 896        |
# | sampleB | locus2 | a      | 584        |
# | sampleB | locus2 | c      | 487        |

gen %>%
  dplyr::select(sample, locus, allele, read_count) %>%
  dplyr::arrange(sample, locus, allele) %>%
  write.csv("output/output2-read-count.csv", row.names = F)

#output3: fasta with all different alleles
#>locus1_a
#ATCT
#>locus1_b
#AACT
#>locus2_a
#CCTC
#>locus2_b
#CCTA

al %>%
  dplyr::mutate(fasta_name = paste(locus, allele, sep = "_")) %>%
  dplyr::arrange(locus, allele) %>%
  {seqinr::write.fasta(
    sequences = as.list(.$seq),
    names = .$fasta_name,
    file.out = "output/output3-alleles.fasta")}

#output4: fasta with alleles per sample
#>sampleA_locus1_a
#ATCT
#>sampleA_locus1_b
#AACT
#>sampleB_locus1_b
#CCTA
#>sampleB_locus1_b
#CCTA
gen %>%
  dplyr::mutate(allele = as.character(allele)) %>%
  dplyr::left_join(al,
                   by = c("allele", "locus")) %>%
  dplyr::mutate(fasta_name = paste(sample, locus, allele, sep = "_")) %>%
  dplyr::arrange(sample, locus, allele) %>%
  {seqinr::write.fasta(
    sequences = as.list(.$seq),
    names = .$fasta_name,
    file.out = "output/output4-alleles-sample.fasta")}
