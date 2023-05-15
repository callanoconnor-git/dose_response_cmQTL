## Run genome scans using qtl2 
## Original Author: Greg Keele
## Adapted by Callan OConnor
## O:12_09_20
## A:10_11_22
options(stringsAsFactors = FALSE)

## Command line arguments
# first argument:    viewer_path         -- path to raw protein/peptide data
# second argument:   dataset             -- dataset list object in QTL Viewer
# third argument:    data_type           -- Form of data to be used: norm or rz
# fourth argument:   output_dir          -- path to write output directory
# fifth argument:    use_kinship         -- logical on whether to use kinship
# sixth argument:    num                 -- number of NAs or 0s
# seventh argument:  key_value           -- NA, zero, or duplicate
# eighth argument:   this_chunk          -- chunk out of 100

## Necessary libraries
library(qtl2)

## Command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if (length(args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    a <- strsplit(args[i], split = '=', fixed = TRUE)[[1]]
    assign(a[1], a[2])
  }
}

## Set number of cores to 4
num_cores <- 4


##Load viewer data
load(viewer_path)
outcome_mat <- readRDS(dataset)

### Check to see if required data are loaded in global environment
stopifnot(c("genoprobs","outcome_mat") %in% ls())

### Pulling pieces from QTL Viewer

## Grab correct outcome_mat column
if (key_value == "NA") {
  counts <- apply(outcome_mat, 2, function(x) sum(is.na(x)))
} else if (key_value == "zero") {
  counts <- apply(outcome_mat, 2, function(x) sum(x == 0))
} else if (key_value == "duplicate") {
  counts <- apply(outcome_mat, 2, function(x) sum(duplicated(x)))
}
this_column <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### Check to see if column exists with the set number of NAs
#stopifnot(!is.na(this_column))
this_column_name<- colnames(as.data.frame(outcome_mat)[this_column])
## Making output file name and stopping if it already exists
output_file <- paste(output_dir,
                     paste0(this_column_name, '.rds'), sep = "/")
stopifnot(!file.exists(output_file))

## Set seed
set.seed(this_chunk)


outputs <- scan1(genoprobs = genoprobs,
                            pheno = outcome_mat[,this_column],
                            addcovar = covar_mat,
                            kinship = K) # 1000 chunks doing 100 permutations (100,000 total)



### Save perm output
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
saveRDS(outputs, output_file)
