
library(here)
library(targets)
library(rhdf5)
library(fs)
library(DropletUtils)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DropletUtils", quietly = TRUE)) BiocManager::install("DropletUtils")


tar_load_globals() # to load paths


dat <- Read10X_h5(proj_paths$cellranger$raw_mtx,use.names = TRUE)


dims <- c(c(sample(1:ncol(dat$`Gene Expression`), size = 10000)),
          c(sample(1:nrow(dat$`Gene Expression`), size = 10000)))

dat_new <- dat$`Gene Expression`[dims]

test_paths <- list()
dataset_name <- "test_01"
test_paths$raw_data_dir <- proj_paths$raw_data_dir
test_paths$cellranger$output_dir <- here(test_paths$raw_data_dir, "cellranger", paste0(dataset_name, "_outs"), "multi")
test_paths$cellbender$output_dir <- here(test_paths$raw_data_dir, "cellbender", dataset_name)
test_paths$vireo$output_dir <- here(test_paths$raw_data_dir, "vireo")

# data file paths
# tag_calls_summary <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","tag_calls_summary.csv"))
# assignment_confidence_table <- read_csv(here(path_cellranger_output_dir,"multiplexing_analysis","assignment_confidence_table.csv"))
test_paths$cellranger$tag_calls_per_cell <- here(test_paths$cellranger$output_dir, "multiplexing_analysis", "tag_calls_per_cell.csv")
test_paths$cellranger$raw_mtx <- here(test_paths$cellranger$output_dir, "count", "raw_feature_bc_matrix.h5")

test_paths$cellbender$output <- here(test_paths$cellbender$output_dir, paste0(dataset_name, "_cellbender_corrected_filtered_seurat.h5"))

#### create parent dirs
paths <- unlist(test_paths, recursive = TRUE, use.names = FALSE)

is_file_like <- (file_exists(paths) & !is_dir(paths)) | nzchar(path_ext(paths))

dirs <- c(path_dir(paths[is_file_like]), paths[!is_file_like])
dirs <- unique(dirs[dirs != "." & !is.na(dirs) & nzchar(dirs)])

dir_create(dirs, recurse = TRUE)


 
# Output file
h5file <- test_paths$cellranger$raw_mtx 

# Create & write groups/datasets
h5createFile(h5file)
#h5createGroup(h5file, "experiment1")
h5write(dat_new, h5file, "Gene Expression")
h5write(dat$`Multiplexing Capture`, h5file, "Multiplexing Capture")

# Always close after writing
H5close()


write

proj_paths$cellbender$output



     
     
here(path_cellbender_output_dir, paste0(dataset_name, "_cellbender_corrected_filtered_seurat.h5"))


proj_paths <-  list()

proj_paths$path_cellranger_output_dir <- 