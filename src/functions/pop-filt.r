#remove monomorphic loci
#input is adegenet genind
#maf is the minor allele frequency threshold. Use "0" for removing monomorphic loci.
#output is filtered genind
remove_monomorphic <- function(genotypes = NULL, maf = 0) {
  if (drop_monomorphic == T) {
  polymorphic_loci <-
    adegenet::locNames(genotypes)[isPoly(genotypes, by = "locus", thres = maf)]
} else if (drop_monomorphic == F) {
  polymorphic_loci <- adegenet::locNames(genotypes)
}
  cat("\nLoci",
  paste(locNames(genotypes)[!(locNames(genotypes) %in% polymorphic_loci)]),
  "were discarded because they were monomorphic\n")
  return(genotypes[loc = polymorphic_loci])
  }

#filter out individuals according to the proportion of alleles genotyped
#input and outputs are genind objects.
filter_ind_missingness <- function(genotypes = NULL) {
    total_alleles <- ploidy * adegenet::nLoc(genotypes)
    alleles_per_sample <- apply(genotypes@tab, 1, sum)
    missingness <- 1 - (alleles_per_sample / total_alleles)
    ind_to_retain <- adegenet::indNames(genotypes)[missingness < max_missing_ind]
    cat("\nSamples",
    paste(indNames(genotypes)[!(indNames(genotypes) %in% ind_to_retain)]),
    "were discarded because they had large amounts of missing data\n")
    return(genotypes[ind_to_retain, ])
}

#filter out loci according to the proportion of alleles genotyped
#input and outputs are genind objects.
missingness_locus <- function(genotypes = NULL) {
  h <-
    tab(genotypes) %>% t() %>%
    as.data.frame() %>%
    dplyr::mutate(locus = locFac(genotypes)) %>%
    group_split(locus, .keep = F) %>%
    plyr::ldply(function(x) colSums(x, na.rm = T) / ploidy) %>%
    rowMeans(na.rm = T) %>% setNames(adegenet::locNames(genotypes))
  loci_to_retain <-
    names(h)[h > (1 - max_missing_loci)]
  cat("\nLoci",
  paste(locNames(genotypes)[!(locNames(genotypes) %in% loci_to_retain)]),
  "were discarded because they had large amounts of missing data\n")
  return(genotypes[loc = loci_to_retain])
  }
#vector with loci non in hw equilibrium

#simple imputation method to impute genotypes to hemizygotes and NA genotypes.
#input and output are genind objects
#warning: it imputs genotypes from randomly resampling genotypes whithin the population.
#the goal of imputation is to overcome limiations at estimating loci not in HWE.
#use at your own responsability. Not reliable for loci with large amount of missing data. It might return innacurate estimates of HWE.

imputation_hw <- function(genind_object) {
  h <-
  pegas::genind2loci(genind_object) %>%
    plyr::ddply(~population, function(x) {
    dplyr::select(x, -population) %>%
      apply(2, function(locusi) {
        locusi[grep(pattern = "/", x = locusi, invert = T)] <- ""
        freq_al <- table(locusi)
        freq_al0 <- freq_al
        freq_al <- freq_al[names(freq_al) != ""]
        allele_vector <- rep(names(freq_al), freq_al)
        no_na <- freq_al0[names(freq_al0) == ""]
        if(length(no_na) > 0) {
          locusi[locusi == ""] <- sample(allele_vector, no_na, replace = T)
        } else {
          locusi
        }
        return(locusi)
      })
  }) %>% as.data.frame() %>%
    dplyr::mutate_all(as.factor) %>%
    pegas::as.loci() %>%
    pegas::loci2genind()
 adegenet::indNames(h) <- adegenet::indNames(genind_object)
 return(h)
}

#remove loci not in HWE. Imputation is necessary. Else, loci with missin calls will be removed.
#inputs and ouputs are genind objects
nonhw <- function(genotypes = NULL, imputate = T) {
  if (imputate == T) genotypes <- imputation_hw(genotypes)
  loci_hwe <-
    pegas::hw.test(genotypes) %>%
    as.data.frame() %>%
    dplyr::filter(Pr.exact > pvalue_hw) %>%
    rownames()
  cat("\nLoci",
  paste(locNames(genotypes)[!(locNames(genotypes) %in% loci_hwe)]),
  "were discarded because they were not in HWE\n")
  return(genotypes[loc = loci_hwe])
}
