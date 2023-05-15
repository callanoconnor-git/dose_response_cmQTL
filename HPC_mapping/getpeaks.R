## Run permutation scans using qtl2 with QTL Viewers
## Original Author: Greg Keele
## Adapted by Callan OConnor
## O:12_09_20
## A:04_04_22
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

library(tidyverse)

load("/pod/2/reinholdt-lab/Callan/perms2/scripts/qtldata.Rdata")

july2022_allpeaks <-
  list.files(path = "/pod/2/reinholdt-lab/Callan/spring_2023/output",
             pattern = "*.rds",
             full.names = T,
             recursive = T) %>% as.data.frame()

library(qtl2)
getpeaks<- function(input){
  a<- readRDS(input)
  qtl_allpeaks <- as.data.frame(find_peaks(scan1_output = a, map = map, threshold=5.0, drop = 1.5)) %>% dplyr::mutate(trait = input)
  list(res=qtl_allpeaks)
}

dfdat <- NULL
for (i in 1:nrow(july2022_allpeaks)) {
  outcome_dats <- getpeaks(input= july2022_allpeaks[i,]) #basespot
  dfdat <- bind_rows(dfdat,
                         outcome_dats$res)
  
}


saveRDS(dfdat, "/pod/2/reinholdt-lab/Callan/spring_2023/output2/qtl_peaksover5.rds")
