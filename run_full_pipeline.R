library(rhdf5)
library(Seurat)
library(ggplot2)
library(tidyr)
library(Matrix)
library(httr)
library(jsonlite)
library(dplyr)
source("get_variants.R") # get_variants_data
source("annotate_variants.R")
source("generate_seurat_object.R")
source("droplet_analysis.R")

FILES<- c("/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_H.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_J.dna+protein.h5",
          #"/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_I.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_A.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_F.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_K.dna+protein.h5",
          "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_G.dna+protein.h5"
          )

seurat_merged <- generate_mega_seurat(FILES)

generate_mega_seurat <- function(FILES){
  seurat_list<- c()
  for(FILE in FILES){
    print(FILE)
    panel_extract<-rhdf5::h5read(FILE, "/assays/dna_read_counts/metadata/panel_name")[1]
    variants_data<- get_variants_data(FILE, min_VAF_cutoff = 1, min_genotyping_rate = 85)
    ggplot(variants_data, aes(x = VAF)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      labs(title = "Distribution of Variant Allele Frequency (VAF)", x = "VAF (%)", y = "Count")
    anns<- annotate_variants(FILE, panel = panel_extract) # get like if coding, exon, and AA change.
    variants_data_to_use <- select_variants(variants_data, anns)
    seurat_obj <- generate_seurat_object(FILE, variants_data_to_use)
    
    seurat_obj <- run_droplet_normalization(seurat_obj, FILE)
    
    sample_name <- unique(seurat_obj@meta.data[["sample_name"]])
    seurat_list[[sample_name]] <- seurat_obj 
  }
  seurat_merged <- merge(seurat_list[[1]], 
                         y = seurat_list[2:length(seurat_list)], 
                         add.cell.ids = names(seurat_list))  # uses list names as prefixes
  seurat_merged <- process_seurat(seurat_merged)
  saveRDS(seurat_merged, "seurat_merged.Rds")
  umap_preview<- DimPlot(seurat_merged, reduction = "umap",label = TRUE, label.size = 6); umap_preview
  umap_preview
  return(seurat_merged)
}


run_droplet_normalization<- function(seurat_obj, FILE){
  droplet_metadata<- get_all_droplet_data(FILE)
  background_droplets<-droplet_metadata%>% dplyr::filter(Droplet_type=="Empty")%>%  dplyr::filter(dna_size<1.5&dna_size>0.15)%>%pull(Cell)
  all_cells_protein_matrix<- h5read(FILE, "all_barcodes/protein_read_counts/layers/read_counts")
  colnames(all_cells_protein_matrix) <- h5read(FILE, "all_barcodes/protein_read_counts/ra/barcode") |> as.character() |> gsub("-1", "", x = _)
  rownames(all_cells_protein_matrix) <- as.character(h5read(FILE, "all_barcodes/protein_read_counts/ca/id"))
  # split into one matrix with empty cells (group_background_droplets) 
  empty_drops_matrix_input<- all_cells_protein_matrix[,background_droplets]
  protein_matrix<- all_cells_protein_matrix[,h5read(FILE, "assays/dna_variants/ra/barcode")]
  isotype <- grep("IgG",rownames(protein_matrix),value=TRUE) # vector of isotype control names
  dsb_norm <- dsb::DSBNormalizeProtein( # remove ambient protien noise reflected in counts from empty droplets
    cell_protein_matrix = protein_matrix, # cell-containing droplet raw protein count matrix
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
  seurat_merged<- FindClusters(object = seurat_merged, # direct graph clustering
                            assay = "DSB_norm",
                            resolution = 0.5,
                            algorithm = 3,
                            verbose = TRUE)
  seurat_merged<- RunUMAP(object = seurat_merged,
                       seed.use = 1333,
                       min.dist = 0.4, 
                       n.neighbors = 30,
                       features=rownames(seurat_merged@assays[["DSB_norm"]]),
                       verbose = TRUE)
  return(seurat_merged)
}




