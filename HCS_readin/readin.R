
setwd('~/Desktop/october/')
options(scipen = 999) # suppress sci notation
# contains harmony plate collection class builder
source("scripts/harmony_utils.R")
# warehouse for functions
source("scripts/custom_tools.R") 
require(dplyr)


selected_dir <- "october"
# Warning - this logic expects a ./temp_files data dir for caching
# Set overwrite to TRUE when underlying data dir has been changed
harmony_collection <- harmony_create_collection(dir = selected_dir,
                                                overwrite = TRUE,
                                                temp_dir = "./temp_files"
                                                #custom_rename = FALSE
)

dim(harmony_collection$all_plates)

saveRDS(harmony_collection,"~/Desktop/october/harmony_collection.rds")  