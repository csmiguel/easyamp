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

#for all loci
sink("output/log.txt", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on src/2-genotyping.R\n")
genotyping <-
  loci %>%
  lapply(function(locus) {
    # determination of ASVs with dada2
    asvs <- asv_determination(locus)
# process asvs
    clean_asvs <-
    # remove spurious ASVs
      clean_lowfrequencyvariants(asvs) %>%
    # rename alleles
      rename_alleles(allele_names = al_names_all) %>%
    # remove paralogues
      remove_nontarget(locus_i = locus) %>%
    # tidy genotypes
      tidy_asvs() %>%
    # remove_extra_alleles
      remove_extra_alleles(locus_i = locus) %>%
    # call homozygous
      call_homozygous(ploidy_i = ploidy) %>%
    # add locus name to data
      add_locus_name(locus_i = locus)
    }) %>%
    setNames(loci_set) %>%
    #remove NULL loci; ie those that have not been genotyped
    {.[!(.) %>% sapply(is.null)]}
sink()
#reformat GENOTYPES
#create data frame with genotypes
genotypes <-
  genotyping %>%
    sapply(function(x){
      x$genotypes
      }) %>%
      do.call(what = cbind) %>%
      {row.names(.) <- all_samples; .} %>%
      as.data.frame() %>%
      setNames(names(genotyping))

tidy_genotypes <-
  genotypes %>% tibble::rownames_to_column("sample") %>%
  reshape2::melt(id = "sample") %>%
  dplyr::rename(locus = variable) %>%
  tidyr::separate(col = value, into = c("allele1", "allele2"), sep  = 1) %>%
  reshape2::melt(id = c("sample", "locus")) %>%
  dplyr::rename(allele = value) %>%
  dplyr::arrange(sample, locus, allele)

# reformat frequency
tidy_freqs <-
  names(genotyping) %>%
  lapply(function(x){
    genotyping[[x]]$frq %>% #data frame with frequencies
      tibble::rownames_to_column(var = "sample") %>%
      reshape2::melt(id = "sample") %>%
      dplyr::rename(allele = variable,
                    count = value) %>%
      dplyr::mutate(locus = x)
  }
  ) %>%
  do.call(what = rbind)

#join genotypes and frequencies
tidy_genotype_frqs <-
  tidy_genotypes %>%
    dplyr::filter(allele != "0") %>%
    dplyr::left_join(tidy_freqs,
                     by = c("sample" = "sample", "locus" = "locus", "allele" = "allele"))

#reformat ALLELE sequences
#create data frame with allele sequences
al_seqs <-
  genotyping %>%
    lapply(function(x) x$allele_seq) %>%
    reshape2::melt() %>%
    dplyr::rename(allele = name,
                  sequence = seq,
                  locus = L1) %>%
    dplyr::select(locus, allele, sequence) %>%
    dplyr::arrange(locus, allele)

#write alleles to fasta
seqinr::write.fasta(sequences = as.list(al_seqs$sequence),
            names = paste(al_seqs$locus, al_seqs$allele, sep = "_"),
            file.out = "output/output3-alleles.fasta")

#write multifasta: all alleles for all individuals
tidy_genotype_frqs %>%
  dplyr::left_join(al_seqs, by = c("locus" = "locus", "allele" = "allele")) %>%
  {seqinr::write.fasta(sequences = as.list(.$sequence),
                    names = paste(.$sample, .$locus, .$allele, sep = "_"),
                    file.out = "output/output4-individuals-x-alleles.fasta")}

#write allele frequencies
write.table(tidy_genotype_frqs,
            file = "output/output2-allele-freqs.txt",
            quote = F,
            row.names = F)

#write genotype table
write.table(genotypes, file = "output/output1-genotypes.tsv", quote = F)

#save objects
saveRDS(genotypes, file = "data/intermediate/genotypes.rds")
saveRDS(al_seqs, file = "data/intermediate/allele_sequences.rds")
saveRDS(genotyping, file = "data/intermediate/genotyping.rds")
