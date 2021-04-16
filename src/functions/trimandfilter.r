# trim and filter reads
# paired-ends
trim_pe <- function(locus = NULL) {
  cat("\nTruncating ", locus, "\n")
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
    NULL
   } else {
  #get optimal truncation lengths
  if (is.null(truncation)) {
  trunc_f <- optimal_truncation_length(fnFs)
  trunc_r <- optimal_truncation_length(fnRs)
} else { # get truncation length from parameters
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
  if (sum(out[, 2]) == 0) cat("\nFiltering removed all reads for", locus, "\n")
  out
  }
}

#single-ends
trim_se <- function(locus = NULL) {
  fnFs <- sort(list.files("data/intermediate",
                          pattern = paste0("\\.", locus, ".1.fastq.gz"),
                          full.names = TRUE))
  sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
  #if there are no files with the locus (ie, not amplified),
  if (length(fnFs) == 0) {
    cat("\nNo files for", locus)
    NULL
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
  if (sum(out[, 2]) == 0) cat("\nFiltering removed all reads for", locus, "\n")
  out
    }
  }

  #determine the most optimal length for the truncation of sequencing reads.
  # The function reads a list of files, computes the length of their sequences
  # and determines which length cut off yields most bases for downstream
  # filtering with dada2::filterAndTrim, using the argument truncLen.
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
