# genotype function: from each locus, it determines asvs, cleans the genotypes and tidies data.
# iterates over errors by using tryCatch, so loci causing problems will output NULL.
genotype <- function(locus) {
  out <- tryCatch(
   expr = {
    # determination of ASVs with dada2
    if (sepe == "pe")
      asvs <- asv_determination_pe(locus)
    if (sepe == "se")
      asvs <- asv_determination_se(locus)
    # process asvs
    clean_asvs <-
      # remove spurious ASVs
      clean_lowfrequencyvariants(asvs1 = asvs) %>%
      # rename alleles
      rename_alleles() %>%
      # remove paralogues
      remove_nontarget(locus_i = locus) %>%
      # tidy genotypes
      tidy_asvs() %>%
      # remove_extra_alleles
      remove_extra_alleles(locus_i = locus) %>%
      # call homozygous
      call_homozygous() %>%
      # add locus name to data
      add_locus_name(locus_i = locus)
      # it will be returned in case there is no condition
      # (e.g. warning or error).
    },
    error = function(cond) {
      message(paste("\nAn error ocurred while genotyping", locus,
      ". This can be caused for the following reasons: no reads passed the quality filters, no data after cleaning genotypes or other errors. Inspect functions whithin 'genotype' and output/log.txt for further details"))
      #if there is any ambigous position
      #if filtering returns no reads, then assign 00 to each genotype"))
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    },
    finally = {
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      message(paste("\nlocus", locus, "processed"))
    }
  )
  return(out)
}

#determination of asvs using dada2:
#importantly:
# setting "pooled = T" in dada2::dada at the same time that having very little
# coverage with bad quality reads can lead to false genotype calls in indivuals
# with rare alleles. "pooled = T" is only recommended if read quality are high
# and coverage is high. Nevertheless, in some ocassion it could yield higher
# sentitivity for detecting rare alleles, if read quality is high.
# Else, it is recommended a more conservative approach and use "pooled = F".
# This option will return, NULL genotype calls for many individuals but more
# reliable genotypes for the positive calls. Parameters in cutadapt (-q)
# and dada2::dada (maxEE) can be tunned to expect higher quality reads.

# locus, locus name as in file with primers
asv_determination_pe <- function(locus = NULL,
                  merge_threshold = get("merge_threshold", envir = .GlobalEnv),
                  truncation = get("truncation", envir = .GlobalEnv)) {
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
   } else {
  #get optimal truncation lengths
  if (is.null(truncation)) {
  trunc_f <- optimal_truncation_length(fnFs)
  trunc_r <- optimal_truncation_length(fnRs)
    } else {
  trunc_f <- truncation[1]
  trunc_r <- truncation[2]
            }
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
  if (sum(out[, 2]) == 0) {
    cat("Filtering removed all reads for", locus)
    } else {

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
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE, pool = F)
  dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE, pool = F)

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
    sapply(function(x) {
      h <- sum(mergers[[x]]$abundance)
      h2 <- sum(dadaFs[[x]]@.Data[[2]]$abundance)
      h / h2
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

#for single end reads
asv_determination_se <- function(locus = NULL,
                        truncation = get("truncation", envir = .GlobalEnv)) {
  cat("\nStarting to genotype", locus, "\n")
  fnFs <- sort(list.files("data/intermediate",
                          pattern = paste0("\\.", locus, ".1.fastq.gz"),
                          full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
  #if there are no files with the locus (ie, not amplified),
  if (length(fnFs) == 0) {
    cat("\nNo files for", locus)
   } else {
  #get optimal truncation lengths
  if (is.null(truncation)) {
  trunc_f <- optimal_truncation_length(fnFs)
    } else {
  trunc_f <- truncation[1]
            }
  # paths to filtered sequences
  filtFs <- file.path(filt_path, paste0(sample.names, "_", locus,
    "_F_filt.fastq.gz"))
    #filter. There cannot be ambiguous positions
  out <- dada2::filterAndTrim(fnFs, filtFs,
                       maxN = 0,
                       maxEE = 3,
                       compress = TRUE,
                       multithread = TRUE,
                       truncLen = trunc_f)
  if (sum(out[, 2]) == 0) {
    cat("Filtering removed all reads for", locus)
    } else {
  #error matrices
  errF <- dada2::learnErrors(filtFs, multithread = TRUE)
  #dereplication
  sample.names.filt <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  derepFs <-
    dada2::derepFastq(filtFs, verbose = TRUE) %>% setNames(sample.names.filt)
  #remove filtered files
  file.remove(filtFs)
  #infer ASV for pooled samples
  dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE, pool = F)
  #Construct sequence table from the provided list of samples:
  seqtab <- dada2::makeSequenceTable(dadaFs) %>%
          dada2::removeBimeraDenovo() #Remove bimeras
          return(seqtab)
        }
      }
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

# helper function of clean_lowfrequencyvariants. It does the filtering.
filter_cov <- function(vector_reads = NULL, cov, ab) {
  h <- max(vector_reads)
  vector_reads[vector_reads / h < ab | vector_reads < cov] <- 0
  vector_reads
  }
clean_lowfrequencyvariants <- function(asvs1 = NULL,
  ab = get("ab", envir = .GlobalEnv),
  cov = get("cov", envir = .GlobalEnv)) {
  #cols are alleles, rows being samples and cells being variant count.
  #removes alleles with read count below min coverage and unbalanced read count.
  #removes alleles with all 0 read_count
  #outputs dataframe with t(input): rows are alleles
  y <- asvs1 %>% t() %>% as.data.frame()
  y[] <- lapply(y, filter_cov, cov, ab)
  #remove alleles with 0 count for all individuals
  y <- y[apply(y, 1, function(x) !all(x == 0)), ] %>% t()
  as.data.frame(y)
  }

#renames alleles according to vector allele_names. It outputs a list with:
# [1] asv table with alleles renamed: individuals in rows and alleles in cols.
# [2] dataframe with different allele names and their sequence.
rename_alleles <- function(asvs1 = NULL,
  al_names_all = get("al_names_all", envir = .GlobalEnv)) {
    al_names <- al_names_all[seq_len(ncol(asvs1))]
    assertthat::assert_that(ncol(asvs1) < length(al_names_all),
      msg = "too many alleles. Consider increasing al_names")
    #create data frame with allele sequences
    al_seq <- data.frame(name = al_names,
                     seq = colnames(asvs1))
    #name alleles
    colnames(asvs1) <- al_names
    list(asvs = asvs1, al_seq = al_seq)
  }

#remove alleles which are suscipious of being non-target sequences
# asvs, is the object generated after rename_alleles
remove_nontarget <- function(asvs = NULL,
  threshold_size = get("threshold_size", envir = .GlobalEnv),
  threshold_clust = get("threshold_clust", envir = .GlobalEnv),
  reference_alleles = get("reference_alleles", envir = .GlobalEnv),
  locus_i = NULL) {
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
    dplyr::select(name) %>%
    dplyr::pull()

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
remove_extra_alleles <- function(asvs1 = NULL,
  ploidy = get("ploidy", envir = .GlobalEnv),
  locus_i = NULL) {
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
  #filter list of alleles
    unique_alleles <- unique(asvs1$asvs$allele)
    asvs1$al_seq <- filter(asvs1$al_seq, name %in% unique_alleles)
    asvs1
    }

#call homozygous GENOTYPES. it calls homozygous genotypes whenever there is only one allele for a diploid locus that overpasses "threshold_homozygosity".
# input is list with tidy asvs and allele sequences.
call_homozygous <- function(asvs1 = NULL,
  ploidy = get("ploidy", envir = .GlobalEnv),
  threshold_homozygosity = get("threshold_homozygosity", envir = .GlobalEnv),
  ...) {
  asvs1$asvs %<>%
  plyr::ddply(~sample, function(x) {
      # create homozygous calls for diploids
      if(ploidy == 2) {
        #call homozygous
        if(nrow(x) == 1 & x$read_count > threshold_homozygosity) {
          x$read_count <- x$read_count / ploidy #if homozygous, divide read_count by ploidy
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
