# Created by use_targets().
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(here)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "reticulate",
    "readr",
    "dplyr",
    "tidyr",
    "purrr",
    "ggplot2",
    "Seurat",
    "tidyseurat",
    "BPCells",
    "stringr",
    "ggridges",
    "forcats",
    "rhdf5",
    "here",
    "patchwork",
    "hdf5r",
    "ggplotify",
    "ggExtra",
    "gridExtra",
    "gt",
    "fs"
  )
)

#### -----


# Set paths for intermediate data




## create directories if not present
if (!file.exists(path_output_plots_root)) {
  dir.create(path_output_plots_root, recursive = T)
  dir.create(path_intermediate_data_root, recursive = T)
  dir.create(path_output_dashboards, recursive = T)
  dir.create(path_output_dashboards, recursive = T)
}

# Directories paths of Raw data
# path_raw_data_file <- here("..", "..", "raw_data", "raw_data_dell", "cubi_cusco")
# path_cellranger_output_dir <- here(path_raw_data_file, "cellranger", paste0(dataset_name, "_outs"), "multi")
# path_cellbender_output_dir <- here(path_raw_data_file, "cellbender", dataset_name)
# path_vireo_output_dir <- here(path_raw_data_file, "vireo")

dataset_name <- "hub_10"
# input directory paths
proj_paths <-  list()

proj_paths$raw_data_dir <- here("..", "..", "raw_data", "raw_data_dell", "cubi_cusco")
proj_paths$cellranger$output_dir <- here(proj_paths$raw_data_dir, "cellranger", paste0(dataset_name, "_outs"), "multi")
proj_paths$cellbender$output_dir <- here(proj_paths$raw_data_dir, "cellbender", dataset_name)
proj_paths$vireo$output_dir <- here(proj_paths$raw_data_dir, "vireo")


# data file paths
# tag_calls_summary <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","tag_calls_summary.csv"))
# assignment_confidence_table <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","assignment_confidence_table.csv"))
proj_paths$cellranger$tag_calls_per_cell <- here(proj_paths$cellranger$output_dir, "multiplexing_analysis", "tag_calls_per_cell.csv")
proj_paths$cellranger$raw_mtx <- here(proj_paths$cellranger$output_dir, "count", "raw_feature_bc_matrix.h5")

proj_paths$cellbender$output <- here(proj_paths$cellbender$output_dir, paste0(dataset_name, "_cellbender_corrected_filtered_seurat.h5"))
proj_paths$cellranger$output <-  here(proj_paths$cellranger$output_dir, )





# change later
# proj_paths$output_root_step <- here("output", dataset_name,step) 
# fs::dir_create()
# proj_paths$plots_root <- here(proj_paths$output_root_step, "plots_and_plot_data")
# proj_paths$dashboards <- here(proj_paths$output_root_step, "dashboards")


#####

leiden_or_louvain <- "louvain"
counts_per_cell_pre_filter_global <- 50

resolution_clustering_1 <- 0.8
resolution_clustering_2 <- 0.2
dim_number_pca_and_neighbors <- 15


hashtag_conversion <- c(
  "CMO301" = "hLOA_BIHi001-B",
  "CMO302" = "hLOA_BIHi005-A",
  "CMO303" = "hLOA_BIHi250-A",
  "CMO304" = "hLOA_UCSFi001-A"
)
#### --------------------------------------


tar_source() # Run the R scripts in the R/

list(
  tar_target(cellbender_output,
    proj_paths$cellbender$output,
    format = "file"
  ),
  tar_target(cellranger_calls,
             proj_paths$cellranger$tag_calls_per_cell,
             format = "file"),
  tar_target(bpcells_cellranger,
    convert_to_bp_cells_object(cellranger_calls, data_slot = "Gene Expression"),
    format = "file"
  ),
  tar_target(bpcells_cellranger_cmo,
    convert_to_bp_cells_object(cellranger_calls, data_slot = "Multiplexing Capture"),
    format = "file"
  ),
  tar_target(bpcells_cellbender,
    convert_to_bp_cells_object(cellbender_output, datasource = "cellbender"),
    format = "file"
  ),
  tar_target(knee_plot_cellranger,
    make_knee(bpcells_cellranger, data_type = "rna"),
    format = "file"
  ),
  tar_target(knee_plot_cellranger_cmo,
    make_knee(bpcells_cellranger_cmo, data_type = "cmo"),
    format = "file"
  ),
  tar_target(
    seurat_obj_raw_joined_cellranger_cellbender,
    create_joined_seurat_obj(
      bpcells_cellbender,
      bpcells_cellranger,
      bpcells_cellranger_cmo,
      counts_per_cell_pre_filter = 50,
      convert_vector = hashtag_conversion
    )
  ),
  tar_target(
    seurat_obj_raw_joined_cellranger_cellbender_vireo,
    add_vireo(
      seurat_obj_raw_joined_cellranger_cellbender,
      path_vireo_output_dir
    )
  ),
  tar_target(setup,
             "0_setup"),
  tar_target(
    seurat_obj_basic_processed,
    basic_seurat_processing(
      obj = seurat_obj_raw_joined_cellranger_cellbender_vireo,
      convert_vector = hashtag_conversion,
      step_no. = setup,
    ) # ,
    # format="file"
  ),
  tar_target(
    seurat_obj_basic_processed_w_metadata,
    add_meta_data(obj = seurat_obj_basic_processed)
  ),
  tar_target(grobs_basic_processed,
    generate_plot_grobs(
      step_no = setup,
      dataset = dataset_name,
      obj = seurat_obj_basic_processed_w_metadata
    ),
    format = "file"
  ),
  tar_target(
    dashboard_1,
    make_dashboard_from_grobs(
      path_grob_filtered = grobs_basic_processed,
      step_no = setup,
      dataset = dataset_name
    )
  )
)
