library(rhdf5)
library(Seurat)
library(ggplot2)
library(tidyr)
library(Matrix)
library(httr)
library(jsonlite)
library(dplyr)
library(HDF5Array)
library(DelayedMatrixStats)
source("get_variants.R") # get_variants_data
source("annotate_variants.R")
source("generate_seurat_object.R")
source("droplet_analysis.R")


FILES<- c("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_I.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_H.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_J.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_A.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_F.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_K.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_G.dna+protein.h5"
          )

seurat_merged <- generate_mega_seurat(FILES, normalize_before_merging= TRUE)


generate_mega_seurat <- function(FILES, normalize_before_merging= TRUE){
  seurat_list<- c()
  for(FILE in FILES){
    print(FILE)
    panel_extract<-rhdf5::h5read(FILE, "/assays/dna_read_counts/metadata/panel_name")[1]
    variants_data<- get_variants_data(FILE, min_VAF_cutoff = 1, min_genotyping_rate = 85)
    anns<- annotate_variants(FILE, panel = panel_extract) # get like if coding, exon, and AA change.
    variants_data_to_use <- select_variants(variants_data, anns)
    seurat_obj <- generate_seurat_object(FILE, variants_data_to_use)
    sample_name <- unique(seurat_obj@meta.data[["sample_name"]])
    #remove poor quality cells 
    seurat_obj <- cell_qc(seurat_obj)
    colnames(seurat_obj) <-paste0(sample_name, "_", colnames(seurat_obj)) # prepend the samplename to the barcodes.
    if (normalize_before_merging == TRUE){
      seurat_obj <- run_droplet_normalization(seurat_obj)
    }
    seurat_list[[sample_name]] <- seurat_obj 
  }
  
  seurat_merged <- merge(seurat_list[[1]], 
                         y = seurat_list[2:length(seurat_list)])
  
  if (normalize_before_merging == FALSE){ # Perform dsb norm on whole cohort
    seurat_merged <- run_droplet_normalization(seurat_merged)
  }
  
  seurat_merged <- process_seurat(seurat_merged) # scale, pca, umap, etc.
  saveRDS(seurat_merged, "seurat_merged.Rds")
  umap_preview<- DimPlot(seurat_merged, reduction = "umap", label = TRUE, label.size = 6); umap_preview
  umap_preview
  return(seurat_merged)
}


# remove cells with average low quality
cell_qc<- function(seurat_obj){
  # Extract the DP matrix
  dp_matrix <- seurat_obj@assays[["NGT"]]@layers[["GQ"]]
  # Calculate a per-cell summary (e.g., column sums or means)
  cell_scores <- colSums(dp_matrix)  # or colMeans()
  # Find the 25th percentile threshold
  q1_threshold <- quantile(cell_scores, 0.25)
  # Identify cells ABOVE the lowest quartile
  cells_to_keep <- cell_scores > q1_threshold
  # Subset the Seurat object
  seurat_obj <- seurat_obj[, cells_to_keep]
  return(seurat_obj)
}

# h5_file_paths can be a list or individual file path
# goal is to use same function for both norm before and after merging seurat objects
run_droplet_normalization<- function(seurat_obj){
  h5_file_paths <- as.list(unique(seurat_obj@meta.data[["file_path"]]))
  empty_drops_matrix_list <- list()
  good_drops_matrix_list  <- list()
  
  for (f in h5_file_paths) {
    sample_name<- unique(seurat_obj@meta.data$sample_name[seurat_obj@meta.data$file_path == f])
    print(paste("getting droplet data for sample", sample_name))
    droplet_metadata <- get_all_droplet_data(f)
    
    background_droplets <- droplet_metadata %>%
      dplyr::filter(Droplet_type == "Empty") %>%
      dplyr::filter(dna_size < 1.5 & dna_size > 0.15) %>%
      pull(Cell)
    
    all_cells_protein_matrix <- h5read(f, "all_barcodes/protein_read_counts/layers/read_counts")
    
    colnames(all_cells_protein_matrix) <- h5read(f, "all_barcodes/protein_read_counts/ra/barcode") |>
      as.character() |>
      gsub("-1", "", x = _)
    
    rownames(all_cells_protein_matrix) <- as.character(
      h5read(f, "all_barcodes/protein_read_counts/ca/id")
    )
    
    # --- Empty drops ---
    empty_drops_matrix <- all_cells_protein_matrix[, background_droplets, drop = FALSE]
    colnames(empty_drops_matrix) <-paste0(sample_name, "_", colnames(empty_drops_matrix))
    
    # Downsample if needed
    set.seed(1333)
    n_cols <- ncol(empty_drops_matrix)
    if (n_cols > 10000) {
      sampled_cols <- sample(seq_len(n_cols), 10000)
      empty_drops_matrix <- empty_drops_matrix[, sampled_cols, drop = FALSE]
    }
    
    # Append as list element 
    empty_drops_matrix_list[[length(empty_drops_matrix_list) + 1]] <- empty_drops_matrix
    
    # good droplets: -------
    protein_matrix <- GetAssayData(seurat_obj, assay="Protein", layer="counts")
    good_drops_matrix_list[[length(good_drops_matrix_list) + 1]] <- protein_matrix
  }
  good_drops_matrix_input <- do.call(cbind, good_drops_matrix_list)
  empty_drops_matrix_input <- do.call(cbind, empty_drops_matrix_list)

  isotype <- grep("IgG",rownames(protein_matrix),value=TRUE) # vector of isotype control names
  dsb_norm <- dsb::DSBNormalizeProtein( # remove ambient protien noise reflected in counts from empty droplets
    cell_protein_matrix = good_drops_matrix_input, # cell-containing droplet raw protein count matrix
    empty_drop_matrix = empty_drops_matrix_input, # empty/background droplet raw protein count
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, #  use isotype controls to define technical components.
    isotype.control.name.vec = isotype # isotype controls (IgG markers) to adjust for non-specific antibody binding
  )
  # remove rows named isotype
  dsb_norm <- dsb_norm[!(rownames(dsb_norm) %in% isotype), ]
  dsb_norm_assay <- CreateAssay5Object(counts = dsb_norm)
  seurat_obj[["DSB_norm"]] <- dsb_norm_assay
  DefaultAssay(seurat_obj) <- "DSB_norm"
  return(seurat_obj)
}

# after merging, do scaling and such
process_seurat<- function(seurat_merged){
  seurat_merged <- Seurat::ScaleData(seurat_merged,assay = "DSB_norm", layer = "counts")
  seurat_merged<- Seurat::RunPCA(seurat_merged, features = rownames(seurat_merged@assays[["DSB_norm"]]))
  seurat_merged<- FindNeighbors(object = seurat_merged, k.param = 30)
  seurat_merged<- FindClusters(object = seurat_merged, assay = "DSB_norm",resolution = 0.5, algorithm = 3,verbose = TRUE)
  seurat_merged<- RunUMAP(object = seurat_merged,seed.use = 1333,min.dist = 0.4,  n.neighbors = 30,features=rownames(seurat_merged@assays[["DSB_norm"]]), verbose = TRUE)
  return(seurat_merged)
}




