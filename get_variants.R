


# we only want to look at varaints that are interesting. those are the not WT ones. 
# WT is specified as VAF lower than like 25 i think. 50 is HET and 100 is HOM. 
get_variants_data<- function(FILE, min_VAF_cutoff = 1,  min_genotyping_rate = 85){
  #ngt <- h5read(FILE, "assays/dna_variants/layers/NGT")  # variants x cells

  ngt <- HDF5Array(FILE, "assays/dna_variants/layers/NGT")  # NOT loaded into RAM
  variant_ids <- h5read(FILE, "assays/dna_variants/ca/id")
  barcodes    <- h5read(FILE, "assays/dna_read_counts/ra/barcode")
  
  dimnames(ngt) <- list(as.character(variant_ids), as.character(barcodes))
  
  variants_data <- data.frame(
    id      = variant_ids,
    WT      = DelayedMatrixStats::rowCounts(ngt, value = 0),
    Het     = DelayedMatrixStats::rowCounts(ngt, value = 1),
    Hom     = DelayedMatrixStats::rowCounts(ngt, value = 2),
    Missing = DelayedMatrixStats::rowCounts(ngt, value = 3),
    ALT     = h5read(FILE, "assays/dna_variants/ca/ALT"),
    CHROM   = h5read(FILE, "assays/dna_variants/ca/CHROM"),
    POS     = h5read(FILE, "assays/dna_variants/ca/POS"),
    QUAL    = h5read(FILE, "assays/dna_variants/ca/QUAL"),
    REF     = h5read(FILE, "assays/dna_variants/ca/REF")
  )
  variants_data <- variants_data %>% 
    dplyr::mutate(VAF = ((Het + Hom * 2) / ((WT + Het + Hom) * 2)) * 100) %>%
    dplyr::mutate(genotyping_rate = ((WT + Het + Hom) / (WT + Het + Hom + Missing)) * 100) %>%
    dplyr::arrange(desc(VAF))
  
  variants_data <- dplyr::filter(variants_data, VAF > min_VAF_cutoff) 
  variants_data <- dplyr::filter(variants_data, genotyping_rate > min_genotyping_rate)
  return(variants_data)
}
