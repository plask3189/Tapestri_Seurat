
select_variants <- function(variants_data, anns){
  variants_data_to_use <- left_join(variants_data, anns[c("CDS_change", "AA_change", "id", "CONSEQUENCE", "Class")], by = "id") %>% filter(!is.na(AA_change)) %>% filter(CONSEQUENCE != "synonymous") 
  variants_data_to_use$AA_change <- make.unique(variants_data_to_use$AA_change)
  rownames(variants_data_to_use) <- variants_data_to_use$AA_change
  return(variants_data_to_use)
}

generate_seurat_object<- function(FILE, variants_data_to_use){
  # Read variant and barcode IDs
  variant_ids <- h5read(FILE, "assays/dna_variants/ca/id")
  barcodes    <- h5read(FILE, "assays/dna_variants/ra/barcode")
  # Subset to variants of interest
  row_idx     <- which(variant_ids %in% variants_data_to_use$id)
  matched_idx <- match(variant_ids[row_idx], variants_data_to_use$id)
  # Read and label NGT matrix
  ngt            <- h5read(FILE, "assays/dna_variants/layers/NGT", index = list(row_idx, NULL))
  rownames(ngt)  <- as.character(variants_data_to_use$AA_change[matched_idx])
  colnames(ngt)  <- as.character(barcodes)
  seurat_obj <- CreateSeuratObject(ngt, assay = "NGT")
  
  seurat_obj@assays[["NGT"]]@meta.data <- as.data.frame(variants_data_to_use[, c(
    "id",
    "AA_change",
    "CDS_change",
    "CONSEQUENCE",
    "Class",
    "VAF",
    "genotyping_rate"
  )])

  # Read additional variant layers
  AF <- h5read(FILE, "assays/dna_variants/layers/AF", index = list(row_idx, NULL))
  DP <- h5read(FILE, "assays/dna_variants/layers/DP", index = list(row_idx, NULL))
  GQ <- h5read(FILE, "assays/dna_variants/layers/GQ", index = list(row_idx, NULL))
  rownames(AF) <-rownames(DP) <-rownames(GQ) <- rownames(ngt)
  colnames(AF) <- colnames(DP) <- colnames(GQ) <- colnames(ngt)
  LayerData(seurat_obj, assay = "NGT", layer = "AF") <- AF
  LayerData(seurat_obj, assay = "NGT", layer = "DP") <- DP
  LayerData(seurat_obj, assay = "NGT", layer = "GQ") <- GQ
  sample_name<- unique(h5read(FILE, "assays/dna_read_counts/ra/sample_name"))
  seurat_obj <- AddMetaData(seurat_obj, sample_name, col.name = "sample_name")
  variants_data_to_use$Gene <- sub("\\..*", "", variants_data_to_use$AA_change)
  seurat_obj <- AddMetaData(seurat_obj, variants_data_to_use$Gene, col.name = "Gene")
  
  # add protein counts ----
  protein_read_counts <- h5read(FILE, "assays/protein_read_counts/layers/read_counts")
  colnames(protein_read_counts) <- as.character(h5read(FILE, "assays/dna_variants/ra/barcode"))
  rownames(protein_read_counts)<-as.character(h5read(FILE, "assays/protein_read_counts/ca/id"))
  protein_assay <- CreateAssay5Object(counts = protein_read_counts)
  seurat_obj[["Protein"]] <- protein_assay
  
  return(seurat_obj)
}




