
#removes fastq files with no reads from a given directory
remove_empty_fastq <- function(path_fastq = NULL) {
  #path_fastq is the directory where the set of fastq files are.
  # foward
  fastq <- sort(list.files(path = path_fastq,
                          pattern = "fastq.gz", full.names = TRUE))
  no_reads <-
    fastq %>%
      sapply(function(x) {
        ShortRead::readFastq(x) %>%
          length()
      }) %>%
      {. == 0}
  file.remove(fastq[no_reads])
  cat(sum(no_reads), "files removed:\n", fastq[no_reads])
}
