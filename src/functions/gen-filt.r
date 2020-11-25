#returns T/F vector for unscored loci (T)
unscored <- function(gen) {
    apply(gen, 2, function(x) {
    sapply(x, function(y) y == "00") %>% all()
    })
  }

#returns T/F vector for monomorphic loci (T)
monomorphic <- function(gen) {
  apply(gen, 2, function(x) {
    #get vector with all different alleles for each locus
    h <- sapply(x, strsplit, split = "") %>%
          unlist() %>% unique() %>% {.[. != "0"]}
    #T/F if monomorphic (only one allele present == length is 1)
    length(h) == 1
    }
    )
  }

#returs T/F vector for individuals having missing data below threshold
ind <- function(gen) {
    apply(gen, 1, function(x) {
    #how many "00"
    h <- grepl("00", x) %>% sum()
    #how many "[a-z]0"
    h0 <- sum(grepl("0", x)) - h
    #proportion of genotype calls with missing data
    p <- (h * 2 + h0) / (length(x) * 2)
    #indivual xx has less missing data than threshold?
    p < max_missing_ind
  })
}

#returs T/F vector for loci having missing data below threshold
locus <- function(gen) {
  apply(gen, 2, function(x) {
    #how many "00"
    h <- grepl("00", x) %>% sum()
    #how many "[a-z]0"
    h0 <- sum(grepl("0", x)) - h
    #proportion of genotype calls with missing data
    p <- (h * 2 + h0) / (length(x) * 2)
    #indivual xx has less missing data than threshold?
    p < max_missing_loci
  })
}

#vector with loci non in hw equilibrium
nonhw <- function(gen = NULL, p = 0.01){
  #genotypes is a data frame with genotypes and pvalue the threshold for the hw test
  gen[, !monomorphic(gen)] %>%
  {adegenet::df2genind(.,
                       ncode = 2,
                       ind.names = row.names(.),
                       loc.names = names(.),
                       sep = "",
                       NA.char = "0",
                       ploidy = 2,
                       type = "codom")} %>%
  pegas::hw.test() %>%
  {rownames(.)[.[, 4] < p]}
}
