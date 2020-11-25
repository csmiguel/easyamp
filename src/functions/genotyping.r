asv_determination <- function(locus = NULL) {
  cat("\nStarting to genotype", locus, "\n")
  #list of fastq files for a given locus
  fnFs <- sort(list.files("data/intermediate",
                          pattern = paste0("\\.", locus, ".1.fastq.gz"),
                          full.names = TRUE))
  fnRs <- sort(list.files("data/intermediate",
                          pattern = paste0("\\.", locus, ".2.fastq.gz"),
                          full.names = TRUE))
#samples from the fastq files
  sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
  #if there are no files with the locus (ie, not amplified),
  if (length(fnRs) == 0 | length(fnFs) == 0) {
    cat("\nNo files for", locus)
  }else {

  #get optimal truncation lengths
  trunc_f <- optimal_truncation_length(fnFs)
  trunc_r <- optimal_truncation_length(fnRs)

  # paths to filtered sequences
  filtFs <- file.path(filt_path, paste0(sample.names, "_", locus,
    "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_", locus,
    "_R_filt.fastq.gz"))

  #filter. There cannot be ambiguous positions
  out <- dada2::filterAndTrim(fnFs, filtFs,
                       fnRs, filtRs,
                       maxN = 0,
                       maxEE = c(3, 6),
                       compress = TRUE,
                       multithread = TRUE,
                       truncLen = c(trunc_f, trunc_r))
  #if reads are too short to overlap
  #if no reads pass the quality criteria
  #if there is any ambigous position
  #if filtering returns no reads, then assign 00 to each genotype
  if (sum(out[, 2]) == 0) {
    cat("Filtering removed all reads for", locus)
    }else {

  #error matrices
  errF <- dada2::learnErrors(filtFs, multithread = TRUE)
  errR <- dada2::learnErrors(filtRs, multithread = TRUE)

  #dereplication
  sample.names.filt <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  derepFs <-
    dada2::derepFastq(filtFs, verbose = TRUE) %>% setNames(sample.names.filt)
  derepRs <-
    dada2::derepFastq(filtRs, verbose = TRUE) %>% setNames(sample.names.filt)
  #remove filtered files
  file.remove(filtFs, filtRs)
  #infer ASV for pooled samples
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE, pool = T)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE, pool = T)

  #Merge paired reads. Spurious sequence variants are further reduced
  #by merging overlapping reads.
  mergers <-
    dada2::mergePairs(dadaFs, derepFs,
                      dadaRs, derepRs,
                      verbose = TRUE, minOverlap = 10, maxMismatch = 0)
  #if reads do not merge, just concanetate them
  # average proportion of reads that have been merged
  prop_merged <-
    seq_along(mergers) %>%
    sapply(function(x){
      h <- sum(mergers[[x]]$abundance)
      h2 <- sum(dadaFs[[x]]@.Data[[2]]$abundance)
      h/h2
    }) %>% mean()
  if (prop_merged < merge_threshold) {
    mergers <-
      dada2::mergePairs(dadaFs, derepFs,
                        dadaRs, derepRs,
                        verbose = TRUE, justConcatenate = T)
  cat(paste("\nReads from", locus, "were concatenated instead of merged."))
  }

  #Construct sequence table from the provided list of samples:
  seqtab <- dada2::makeSequenceTable(mergers) %>%
          dada2::removeBimeraDenovo() #Remove bimeras
          return(seqtab)
        }
      }
    }

clean_lowfrequencyvariants <- function(asvs1 = NULL, ab1 = ab, cov1 = cov) {
  #table_ab is a matrix, array or dataframe with cols being alleles,
  #rows being samples and cells being variant count.
  #removes alleles with read count below min coverage and unbalanced read count.
  apply(asvs1, 1, function(x) {
    h <- max(x)
    #set to 0 the count of alleles with a coverage below threshold
    #set to 0 unbalance alleles: read count below a minor allele frequency threshold.
    y <- x
    y[x / h < ab1 | x < cov1] <- 0
    y
  }) %>%
  #remove alleles with 0 count for all individuals
  {.[apply(., 1, sum) > 0, ]} %>%
    t() %>%
    as.data.frame()
  }

# determine the most optimal length for the truncation of sequencing reads. The function reads a list of files, computes the length of their sequences and determines which length cut off yields most bases for downstream filtering with dada2::filterAndTrim, using the argument truncLen.
optimal_truncation_length <- function(fq_files) {
  fq_files %>%
    lapply(function(x){
    ape::read.fastq(x) %>%
      sapply(length)
    }) %>%
    unlist() %>% table() %>% as.data.frame() %>%
    dplyr::rename(len = 1) %>%
    dplyr::mutate(len = as.numeric(as.character(len))) %>%
    dplyr::arrange(desc(len)) %>%
    dplyr::mutate(cum_freq = cumsum(Freq),
                  cum_nt = cum_freq * len) %>%
    dplyr::arrange(desc(cum_nt)) %>%
    dplyr::pull(len) %>% {.[1]}
  }

#renames alleles according to vector allele_names. It outputs a list with:
# [1] asv table with alleles renamed: individuals in rows and alleles in cols.
# [2] dataframe with different allele names and their sequence.
rename_alleles <- function(asvs1 = NULL, allele_names = NULL) {
    al_names <- allele_names[seq_len(ncol(asvs1))]
    assertthat::assert_that(ncol(asvs1) < length(allele_names),
      msg = "too many alleles. Consider increasing al_names")
    #create data frame with allele sequences
    al_seq <- data.frame(name = al_names,
                     seq = names(asvs1))
    #name alleles
    colnames(asvs1) <- al_names
    list(asvs = asvs1, al_seq = al_seq)
  }

#remove alleles which are suscipious of being non-target sequences
# asvs, is the object generated after rename_alleles
# reference_alleles, is the path to a fasta file with known sequences of real alleles. The name of the correspondng fasta sequence must be identical to the locus name in the object "loci".
# if reference_alleles are provided, the corresponding locus is taken a as reference for determining the target sequences in "asvs". Otherwise, the allele more frequent across individuals is taken as a reference (as being true target sequence).
remove_nontarget <- function(asvs = NULL,
  reference_alleles = NULL, locus_i = locus) {
  #set the reference allele.
  if (!is.null(reference_alleles)) {
    refs <- ape::read.FASTA(reference_alleles)
    names(refs) <- tolower(names(refs))
    homol_al <- refs[locus_i][1]
    cat(paste("\nA reference allele for", locus_i,
     " was provided as a reference for removing non target alleles\n"))
  } else {
  #the most commong allele is considered the homologue
  homol_al <-
    apply(asvs$asvs, 2, function(x) sum(x > 0)) %>%
    sort(decreasing = T) %>%
    names() %>% .[1] %>%
    {dplyr::filter(asvs$al_seq, name == .)} %>%
    dplyr::pull(2) %>%
    {lapply(strsplit(., ""), tolower)} %>%
    setNames(locus_i) %>%
    ape::as.DNAbin()
    cat(paste("\nThe most frequent allele for", locus_i,
    "was taken as the reference allele for removing non target alleles\n"))
  }

  # alleles not differing too much in size (threshold_size)
  # length of homologous allele
  homol_length <- length(homol_al[[1]])
  al_kept_size <-
    asvs$al_seq %>%
    dplyr::mutate(seq_length = nchar(seq)) %>%
    dplyr::mutate(prop_length =
                    abs((seq_length - homol_length) / homol_length)) %>%
    dplyr::filter(prop_length < threshold_size) %>%
    dplyr::select(name) %>% dplyr::pull()

  # alleles with sequences not differing too much from homologue (threshold_clust).
  # the sensitiviy can be tunned changing "threshold_clust", or "k" argument from kdistance.
  # some values for guidance:
  # k=6 and threshold_clust = 0.12, corresponds to p-distances of around 10%.
  # k=6 and threshold_clust = 0.23, corresponds to p-distances of around 20%.
  # k=6 and threshold_clust = 0.65, correspond to random DNA sequences.
  al_kept_seq <-
      lapply(strsplit(asvs$al_seq$seq, ""), tolower) %>%
      setNames(asvs$al_seq$name) %>%
      ape::as.DNAbin() %>%
      c(homol_al) %>%
      kmer::kdistance(k = 6) %>%
      as.matrix() %>% as.data.frame() %>%
      dplyr::select(eval(locus_i)) %>%
      tibble::rownames_to_column("allele") %>%
      dplyr::filter(allele != eval(locus_i)) %>%
      dplyr::rename(homol = 2) %>%
      dplyr::filter(homol < threshold_clust) %>%
      dplyr::select(allele) %>% dplyr::pull()
  #list of alleles to keep
  kept_al <- unique(c(al_kept_size, al_kept_seq))
  al_droped <- asvs$al_seq$name[!(asvs$al_seq$name %in% kept_al)]
    if(length(al_droped) == 0) al_droped <- "No"
  #remove alleles from genotype matrix
  asvs$asvs %<>% dplyr::select(all_of(kept_al))
  asvs$al_seq %<>% dplyr::filter(name %in% kept_al)
  cat(al_droped, "alleles from", locus_i,
  "were/was dropped because of being suspicious of not being target sequences\n")
  return(asvs)
    }

#tidy asvs in tidy format
tidy_asvs <- function(asvs1 = NULL) {
  asvs1$asvs <-
    reshape2::melt(as.matrix(asvs1$asvs)) %>%
      setNames(c("sample", "allele", "read_count")) %>%
      dplyr::filter(read_count > 0)
  asvs1
  }

# flat asvs: remove alleles above if they overpass the expected ploidy. Some loci for a given individual might have more alleles than expected by their ploidy. In order not to interrupt the workflow, it is assumed the alleles with deeper coverage are the right ones. All others are dropped and assumed not to be target sequences. A log is emmitted.
# input is list with tidy asvs and allele sequences.
remove_extra_alleles <- function(asvs1 = NULL, locus_i = NULL) {
  asvs1$asvs %<>%
  plyr::ddply(~sample, function(x) {
    #1. remove alleles to fit ploidy
    if (nrow(x) > ploidy) {
      x <- dplyr::arrange(x, desc(read_count)) %>%
        {al_droppedx <<- .$allele[3:nrow(x)]; .} %>%
        {.[1:ploidy, ]}
    cat("\n",
        paste(locus_i, "from", unique(x$sample),
        "has more alleles for some samples than the expected ploidy. This can be caused by the co-amplification of other loci (eg. paralogues). Allele(s):\n",
        paste(filter(asvs1$al_seq, name == al_droppedx)$seq, sep = " and "),
        "\nhave been dropped because they had less coverage.\n"))
    } else {x}
    return(x)
  }) #end ddply
    unique_alleles <- unique(asvs1$asvs$allele)
    asvs1$al_seq <- filter(asvs1$al_seq, name %in% unique_alleles)
    asvs1
    }

#call homozygous GENOTYPES. it calls homozygous genotypes whenever there is only one allele for a diploid locus that overpasses "threshold_homozygosity".
# input is list with tidy asvs and allele sequences.
call_homozygous <- function(asvs1 = NULL, ploidy_i = NULL) {
  asvs1$asvs %<>%
  plyr::ddply(~sample, function(x) {
      # create homozygous calls for diploids
      if(ploidy_i == 2) {
        #call homozygous
        if(nrow(x) == 1 & x$read_count > threshold_homozygosity) {
          x$read_count <- x$read_count / ploidy_i #if homozygous, divide read_count by ploidy
          rbind(x, x)
          } else {x}
        } else {
          x
        }
    }) #end of ddply
    asvs1
  }

# add locus name to data in list
add_locus_name <- function(asvs1 = NULL, locus_i = NULL) {
  asvs1$asvs$locus <- locus_i
  asvs1$al_seq$locus <- locus_i
  asvs1
}
