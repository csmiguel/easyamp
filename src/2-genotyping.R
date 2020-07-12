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

#load parameters for genotype calling
source("src/parameters/parameters.r")
#load function to cleas spurious ASVs
source("src/functions/remove-spurious-ASVs.r")

#samples
all_samples <- readLines("data/intermediate/samples-list")
#loci
loci <- readLines("data/intermediate/loci")
#if any subset is desired
loci_set <- loci
#for all loci
genotyping <-
  loci_set %>%
  lapply(function(locus) {
    fnFs <- sort(list.files("data/intermediate",
                            pattern = paste0("*", locus, ".*1.fastq.gz"),
                            full.names = TRUE))
    fnRs <- sort(list.files("data/intermediate",
                            pattern = paste0("*", locus, ".*2.fastq.gz"),
                            full.names = TRUE))
  #samples from the fastq files
    sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
    #if there are no files with the locus (ie, not amplified),
    if (length(fnRs) == 0 | length(fnFs) == 0) {
      cat("no files for", locus)
    }else{
    #create path for sequences that will be filtered downstream
    filt_path <- "data/intermediate/filtered"
    if (!file_test("-d", filt_path)) dir.create(filt_path)

    filtFs <- file.path(filt_path, paste0(sample.names, locus,
      "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, locus,
      "_R_filt.fastq.gz"))

    #filter. There cannot be ambiguous positions
    out <- dada2::filterAndTrim(fnFs, filtFs,
                         fnRs, filtRs,
                         maxN = 0, maxEE = c(3, 6),
                         compress = TRUE, multithread = TRUE)
    #if filtering returns no reads, then assign 00 to each genotype
    if (sum(out[, 2]) == 0) {
      cat("filtering removed all reads for", locus)
      }else{
    #error matrices
    errF <- dada2::learnErrors(filtFs, multithread = TRUE)
    errR <- dada2::learnErrors(filtRs, multithread = TRUE)

    #dereplication
    derepFs <-
      dada2::derepFastq(filtFs, verbose = TRUE) %>% setNames(sample.names)
    derepRs <-
      dada2::derepFastq(filtRs, verbose = TRUE) %>% setNames(sample.names)

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
    if (mergers %>% sapply(function(x) nrow(x) == 0) %>% all()){
      mergers <-
        dada2::mergePairs(dadaFs, derepFs,
                          dadaRs, derepRs,
                          verbose = TRUE, justConcatenate = T)
    }

    #Construct sequence table from the provided list of samples:
    seqtab <- dada2::makeSequenceTable(mergers) %>%
            dada2::removeBimeraDenovo() #Remove bimeras
    #remove spurious ASVs
    hh <- clean_lowfrequencyvariants(seqtab)
    #allele names
    al_names <- letters[1:ncol(hh)]
    #create data frame with allele sequences
    al_seq <- data.frame(name = al_names,
                     seq = names(hh))
    #name alleles
    colnames(hh) <- al_names
    #Determine home/hetero/hemi- zygotes from allele frequency. Output is 2 column data frame.
    geno <-
      hh %>%
        apply(1, function(x) {
          y <- names(x)[x > 0]
          if (length(y) == 0 | length(y) > 2) "00"
          else if (length(y) == 1) {
            if((max(x) - cov) < 5) paste0(y, "0", collapse = "") else paste0(y, y, collapse = "") #for homozygous
          } else paste0(y, collapse = "") #for heterozygous
        }) %>%
        #I will convert it to a data frame and do a left join with all samples
        #in case there is anyone missing
        as.data.frame() %>%
        tibble::rownames_to_column(var = "id") %>%
        {dplyr::left_join(x = data.frame(id = all_samples), y = ., by = "id")} %>%
        tibble::column_to_rownames(var = "id")
    geno[is.na(geno), 1] <- "00"
    list(genotypes = geno, allele_seq = al_seq, frq = hh)
    }
  }
    }) %>%
    setNames(loci_set) %>%
    #remove NULL loci; ie those that have not been genotyped
    {.[!(.) %>% sapply(is.null)]}

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
