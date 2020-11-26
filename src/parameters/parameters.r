#input data
sepe = "pe" #"se" or "pe"; paired-end or single-end reads
#genotype calls
truncation = c(270, 200) # if NULL, the optimal_truncation_length function is run to find the truncation length which leads to the retention of most nucleotides. However, in many cases users will prefer to set the truncation length manually with a vector (eg c(200, 180)) with the desired truncation lengths from F and R reads, in order to ensure sufficient overlap for merging F and R reads.
ab = 0.2 #threshold for allele balance (eg. 0.3 implies that reads_min_allele/reads_max_allele > 0.3)
cov = 3 #minimum coverage for an allele in an individual not be be dropped.
merge_threshold = 0.1 #proportion of reads merged below which concatenation is done.
ploidy = 2 #markers ploidy: either 1 or 2
threshold_homozygosity = 5 #read_count above this threshold is required to call an individual hemozygous. Else, it will be hemizygous.
threshold_size = 0.3 #if (size_ref_allele-size_observed_allele)/size_ref_allele > threshold_size, the observed allele is discarded.
threshold_clust = 0.1 # threshold of kdistances from kmer of a reference target allele to a genotyped allele. Alleles more disimilar than this threshold are considered non-target and discarded.
reference_alleles = NULL # reference_alleles, is the path to a fasta file with known sequences of real alleles. The name of the correspondng fasta sequence must be identical to the locus name in the object "loci". if reference_alleles are provided, the corresponding locus is taken a as reference for determining the target sequences in "asvs". Otherwise, the allele more frequent across individuals is taken as a reference (as being true target sequence).

#pop gen filtering
max_missing_ind = 0.3 #max. proportion of missing calls per individual (recommended 0.3)
max_missing_loci = 0.3 #max. proportion of missing calls per locus (recommended 0.3)
pvalue_hw = 0.05 #threhold pvalue in pegas::hw.test to retain loci
ind_to_remove = NULL # vector with samples to remove from the genotypes in the filtering. For instance, for population genetics data you will have to remove samples from other species in the genotypes.
