clean_lowfrequencyvariants <- function(table_ab) {
  #table_ab is a matrix, array or dataframe with cols being alleles,
  #rows being samples and cells being variant count.
  h <-
    apply(table_ab, 1, function(x){
    h <- max(x)
    x[x/h < ab | x < cov] <- 0
    x
  }) %>% as.data.frame() %>% setNames(colnames(table_ab))
  #when there is only one allele
  if(ncol(table_ab) > 1) h <- t(h)
  #after initial cleaning, some alleles are spurious, and must be removed
  #for that reaseon. I will remove cols which are all 0
  h[, apply(h, 2, sum) > 0] %>%
    as.data.frame(row.names = rownames(table_ab)) %>%
      {setNames(object = ., colnames(table_ab)[1:ncol(.)])}
}
