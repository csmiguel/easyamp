#pop gen filtering
max_missing_ind = 0.3 #max. proportion of non genotyped alleles per individual (recommended 0.3). If 1, no filter is applied.
max_missing_loci = 0.3 #max. proportion of non genotyped alleles per locus (recommended 0.3). If 1, no filter is applied.
pvalue_hw = 0.005 #threhold pvalue in pegas::hw.test to retain loci. If 0, no filter is applied.
drop_monomorphic = T #if false, remove_monomorphic is not applied
