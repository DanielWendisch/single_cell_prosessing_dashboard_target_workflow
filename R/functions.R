# quarto_render_move <- function( #https://github.com/jhelvy/jph/blob/master/R/quarto_render_move.R
#     input,
#     output_file = NULL,
#     output_dir = NULL,
#     ...
# ) {
#   
#   # Get all the input / output file names and paths
#   x <- quarto::quarto_inspect(input)
#   output_format <- names(x$formats)
#   output <- x$formats[[output_format]]$pandoc$`output-file`
#   if (is.null(output_file)) { output_file <- output }
#   input_dir <- dirname(input)
#   if (is.null(output_dir)) { output_dir <- input_dir }
#   output_path_from <- file.path(input_dir, output)
#   output_path_to <- file.path(output_dir, output_file)
#   
#   # Render qmd file to input_dir
#   quarto::quarto_render(input = input, ... = ...)
#   
#   # If output_dir is different from input_dir, copy the rendered output
#   # there and delete the original file
#   if (input_dir != output_dir) {
#     
#     # Try to make the folder if it doesn't yet exist
#     if (!dir.exists(output_dir)) { dir.create(output_dir) }
#     
#     # Now move the output to the output_dir and remove the original output
#     file.copy(
#       from = output_path_from,
#       to = output_path_to,
#       overwrite = TRUE
#     )
#     file.remove(output_path_from)
#     
#     # If the output_dir is the same as input_dir, but the output_file
#     # has a different name from the input file, then just rename it
#   } else if (output_file != output) {
#     file.rename(from = output_path_from, to = output_path_to)
#   }
# }
# 




################################
make_path_seurat_object <- function(step.,
                                    dataset.,
                                    path_interdat_root.){
  
  path_seurat_obj<- here(path_interdat_root.,
                         paste("seurat_obj",
                               step.,
                               dataset.,
                               "rds",
                               sep = "."))
  return(path_seurat_obj)
  
}
####################

convert_to_bp_cells_object <- function(path_output,
                                                  data_slot="Gene Expression",
                                                  datasource="cellranger"){
  dat <- Read10X_h5(path_output,use.names = TRUE)
  if (data_slot=="Gene Expression") { #todo change names to janitor::make_names()
    # todo change to  case_when() syntax
    path_bp_cells <- here("BPcell_matrices", paste0("cellranger_", dataset_name))
  }
  if (data_slot=="Multiplexing Capture") {
    path_bp_cells <- here("BPcell_matrices", paste0("cellranger_cmo_", dataset_name))
  }
  if (datasource=="cellbender") {
    path_bp_cells <- here("BPcell_matrices", paste0("cellbender_", dataset_name))
    if (is.null(names(dat))) {write_matrix_dir(dat,path_bp_cells, overwrite = T)
    } else{
      write_matrix_dir(dat$`Gene Expression`, path_bp_cells, overwrite = T)
    }
    return(path_bp_cells)
    
  }else{write_matrix_dir(mat = dat[[data_slot]],
                   dir = path_bp_cells,
                   overwrite = T)
  return(path_bp_cells)}
}


make_knee <- function(path_bp_cells, data_type="rna"){
  dat <- open_matrix_dir(dir = path_bp_cells)  
  testingthis <- 2
  rna_reads_per_cell <- dat |> colSums()
  knee_plot_rna <- plot_read_count_knee(rna_reads_per_cell, cutoff = 2e3)+
    plot_read_count_knee(rna_reads_per_cell, cutoff = 1e3)+
    plot_read_count_knee(rna_reads_per_cell, cutoff = 5e2)+
    plot_read_count_knee(rna_reads_per_cell, cutoff = 1e2)+
    plot_annotation(title="RNA reads per cell")
  ggsave(knee_plot_rna,
         filename =   here("output",
                           dataset_name,
                           "qc_process",
                           "plots_and_plot_data",
                           paste(dataset_name ,data_type,"read_count_knee.png", sep = "_"))
 )
  
  
}




create_joined_seurat_obj <- function(
    path_bp_cells_cellbender,
    path_bp_cells_cellranger,
    path_bp_cells_cellranger_cmo,
    counts_per_cell_pre_filter=counts_per_cell_pre_filter_global,
    convert_vector
    ){

  dat_ranger <- open_matrix_dir(dir = path_bp_cells_cellranger)
  rna_reads_per_cell <- dat_ranger |> colSums()
  
  dat_cmo <- open_matrix_dir(dir = path_bp_cells_cellranger_cmo)
  dat_ranger <- dat_ranger[,(rna_reads_per_cell>counts_per_cell_pre_filter)]
  dat_cellbender <- open_matrix_dir(dir = path_bp_cells_cellbender)
  
  
  ############ #inital trimming of data
  
  #cells that are present in celranger counts &cellranger hashtag & cellbender output
  intersecting_cells <- intersect(colnames(dat_cellbender), colnames(dat_ranger))
  intersecting_cells <- intersect(intersecting_cells, colnames(dat_cmo))
  
  
  seurat_obj <- CreateSeuratObject(dat_ranger[, intersecting_cells])
  assay_cellbender <- CreateAssay5Object(dat_cellbender[, intersecting_cells])
  
  seurat_obj[["cellbender_RNA"]] <- assay_cellbender
  
  dat_cmo <- dat_cmo[, intersecting_cells]
  dat_cmo_tbl <- as.data.frame(t(as.matrix(dat_cmo)))
  dat_cmo <- CreateAssay5Object(dat_cmo)
  
  seurat_obj[["hashtag_oligos"]] <- dat_cmo
  
  # add hashtag counts to metadata
  dat_cmo_tbl <- dat_cmo_tbl[,names(convert_vector)]
  colnames(dat_cmo_tbl) <- paste("hashtag_id",convert_vector, sep="_")
  seurat_obj <- AddMetaData(object = seurat_obj, dat_cmo_tbl)

  return(seurat_obj)
  }



add_vireo <- function(obj, path_vireo, dataset=dataset_name){
  path_vireo <- here(path_vireo, paste0(dataset, "_vireo_donor_ids.tsv"))
  if (file.exists(path_vireo)) {
    vireo_ids <- read_tsv(path_vireo)
    
    vireo_ids <- vireo_ids |> select(cell, donor_id) |>
      right_join(tibble(cell=colnames(obj))) |>
      mutate(donor_id=if_else(is.na(donor_id),
                              "no vireo data",
                              donor_id))
    
    vireo_ids_vec <- pull(vireo_ids, donor_id)
    names(vireo_ids_vec) <- pull(vireo_ids, cell)
    obj <- AddMetaData(obj,
                       col.name ="vireo_snp_demux",
                       vireo_ids_vec)
    
    obj <- obj |> mutate(vireo=case_when(vireo_snp_demux == "unassigned" ~ "vireo.unassigned" ,
                                         is.na(vireo_snp_demux) ~ "vireo.no data",
                                         TRUE ~ "vireo.cell"))
    
  }else { # some weird error when using else
    obj$vireo <- "not_computed"
    obj$vireo_snp_demux <- "not_computed"}
    return(obj)
    }



basic_seurat_processing <- function(obj,
                                    convert_vector,
                                    dataset=dataset_name,
                                    step_no. ){
  browser()
  for (assay_iter in Assays(obj)) {
    print("new assay #######################################################################")
    print("          #######################################################################")
    print(assay_iter)

    if (assay_iter=="hashtag_oligos") {dimensions <- 1:length(convert_vector)}else{dimensions <- c(1:dim_number_pca_and_neighbors)}

    DefaultAssay(obj) <- assay_iter
    obj <- obj |>
      NormalizeData()

    print("scale data")
    obj <- obj |>
      FindVariableFeatures(assay = assay_iter) |>
      ScaleData()

    print("runPCA")
    obj[[assay_iter]]$scale.data  <- obj[[assay_iter]]$scale.data %>% write_matrix_memory(compress=FALSE)
    obj <- obj |>
      RunPCA(assay = assay_iter, reduction.name = paste0("pca.", assay_iter),npcs =max(dimensions))

    print("neighbors")
    obj <- obj |> FindNeighbors(
      dims = dimensions,
      reduction = paste0("pca.", assay_iter),
      assay = assay_iter)
    gc()

    print("clustering no 1")
    obj <- obj |> FindClusters(
      resolution = resolution_clustering_1,
      verbose = TRUE,
      cluster.name = paste(assay_iter, "clusters", leiden_or_louvain, "res." , resolution_clustering_1, sep = "_"),
      algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
      # graph.name = paste0("kNN_","pca.",assay_iter)
    )

    print("clustering no 2")
    obj <- obj |> FindClusters(
      resolution = resolution_clustering_2,
      verbose = TRUE,
      cluster.name = paste(assay_iter, "clusters", leiden_or_louvain, "res." , resolution_clustering_2, sep = "_"),
      algorithm = ifelse(leiden_or_louvain=="louvain",1,4)
      # graph.name = paste0("kNN_","pca.",assay_iter)
    )


    obj <- obj |>
      RunUMAP(
        dims = dimensions,
        reduction = paste0("pca.", assay_iter),
        assay = assay_iter,
        reduction.name = paste0("umap.", assay_iter)
      )

    gc()
  }
  
  DefaultAssay(obj) <- "RNA"
  path_intermediate_data_root <- here("intermediate_data", dataset_name, step_no. )
  
  if (!dir.exists(path_intermediate_data_root)) {
    dir.create(path_intermediate_data_root, recursive = TRUE)
  }
  rds_path <- make_path_seurat_object(step.=step_no.,
                                      dataset.=dataset,
                                      path_interdat_root.=path_intermediate_data_root)
  # obj |>  write_rds(rds_path)
  # 
  # return(rds_path)
  return(obj)
}








add_meta_data <- function(tags=tag_calls_per_cell_path,
                          obj,
                          counts_per_cell_pre_filter=counts_per_cell_pre_filter_global,
                          convert_vector=hashtag_conversion
                          ){
  
  
  cellbender_probs <- here(path_cellbender_output_dir, paste0(dataset_name,"_cellbender_corrected.h5"))
  tags<- read_csv(tag_calls_per_cell_path)
  # make metadata column: cell_line_cellranger
  #-- match CMOxx to cell line name, adjusted feature_call column from cellranger multi
  cell_barcodes <- tags |>
    right_join(tibble(cell_barcode = colnames(obj)), keep = T) |>
    mutate(cell_line_cellranger = str_replace_all(feature_call, convert_vector) |>
             replace_na("no cell"))
  
  
  # METADATA: singlet_doublet_cellranger
  #-- cellranger calls with levels doublet, singlet no cell(not called a cell from cellranger and therefore not in file per_sample_outs but in file multi>count)
  # adjusted feature_call column from cellranger multi
  
  cell_barcodes <- cell_barcodes |> mutate(singlet_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet",
                                                                               ifelse(cell_line_cellranger == "no cell", "no cell", "singlet")))
  # METADATA: cell_line_doublet_cellranger
  #--- like singlet_doublet_cellranger but with cell line names
  cell_barcodes <- cell_barcodes |> mutate(cell_line_doublet_cellranger = ifelse(grepl("\\|", cell_line_cellranger), "doublet", cell_line_cellranger))
  
  # add relevant dataframe columns to Seurat Metadata
  cell_barcodes <- cell_barcodes |>
    as.data.frame()
  rownames(cell_barcodes) <- cell_barcodes$cell_barcode.y
  cell_barcodes <- cell_barcodes[, c("cell_line_cellranger", "singlet_doublet_cellranger", "cell_line_doublet_cellranger")]
  obj <- AddMetaData(obj, cell_barcodes)
  
  
  
  ### add cellbender cell probabilities ----
  #--------------------------------------------------------
  
  corrected <- H5Fopen(cellbender_probs)
  cell_prob_tbl <- tibble(
    cell = corrected$metadata$barcodes_analyzed,
    cell_probability = corrected$droplet_latents$cell_probability)
  
  H5Fclose(corrected)
  
  cell_prob_tbl <- cell_prob_tbl |> filter(cell %in% colnames(obj))
  cell_prob_vec <- cell_prob_tbl |> pull(cell_probability)
  names(cell_prob_vec) <- cell_prob_tbl |> pull(cell)
  
  # METADATA cellbender_prob_to_be_cell
  # from cellbender output  H5Fopen cell = corrected$metadata$barcodes_analyzed, cell_probability = corrected$droplet_latents$cell_probability
  obj <- AddMetaData(obj, col.name = "cellbender_prob_to_be_cell", cell_prob_vec)
  
  # METADATA cellbender_p95
  # METADATA vireo
  # METADATA cellranger
  
  obj <- obj |> 
    mutate(cellbender_p95=if_else(cellbender_prob_to_be_cell > 0.95, "cellbender.higher0.95" , "cellbender.lower0.95")) |> 
    
    mutate(cellranger=paste("cellranger", singlet_doublet_cellranger, sep="."))
  
  # METADATA cellranger percent_mito
  obj[["percent_mito"]] <- PercentageFeatureSet(object = obj, pattern = "^MT-", assay = "RNA")
  
  
  
  obj@misc$filtering_steps <- list()
  obj@misc$filtering_steps$qc_setup <-paste("counts_per_cell_pre_filter >" ,as.character(counts_per_cell_pre_filter))
  
  
  obj$organoid <-  obj |> separate(cell_line_cellranger, sep = "_", into = "organoid")|> pull(organoid)
  
  obj$cell_line <-  obj |> separate(cell_line_cellranger, sep = "_", into = c("organoid", "cell_line")) |> pull(cell_line)
  return(obj)
  
  
}














