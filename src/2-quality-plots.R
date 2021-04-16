###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: fastq quality plots
#PROJECT: amplicon-genotyping
###.............................................................................
library(ggplot2)
library(dada2)
library(dplyr)
library(ShortRead)

source("src/functions/remove_empty_fastq.r")

remove_empty_fastq(path_fastq = "data/intermediate")
#fastq file paths
# foward
fnFs <- sort(list.files(path = "data/intermediate",
                        pattern = "1.fastq.gz", full.names = TRUE))
# reverse
fnRs <- sort(list.files(path = "data/intermediate",
                        pattern = "2.fastq.gz", full.names = TRUE))

#quality plots
# forward
qualplt_f <- dada2::plotQualityProfile(sample(fnFs, size = 30, replace = T),
                                       n = 100, aggregate = T)
# reverse
qualplt_r <- dada2::plotQualityProfile(sample(fnRs, size = 30, replace = T),
                                       n = 100, aggregate = T)
#add titles
qualplt_f <-
  qualplt_f +
  ggtitle("Forward reads")
qualplt_r <-
  qualplt_r +
  ggtitle("Reverse reads")
#write to file
g <- gridExtra::arrangeGrob(qualplt_f, qualplt_r, nrow = 2)
ggplot2::ggsave("output/quality_plot.pdf", plot = g, width = 4, height = 5)
