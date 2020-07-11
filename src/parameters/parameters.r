#input data
sepe = NULL #"se" or "pe"; paired-end or single-end reads#not yet functional
#genotype calls
ab = 0.2 #threshold for allele balance (eg. 0.3 implies that reads_min_allele/reads_max_allele > 0.3)
cov = 5 #minimum coverage for an allele in an individual not be be dropped.
#pop gen filtering
remove_monomorphic = "yes" #remove remove monomorphic loci
max_missing_ind = 0.3 #max. proportion of missing calls per individual (recommended 0.3)
max_missing_loci = 0.3 #max. proportion of missing calls per locus (recommended 0.3)
