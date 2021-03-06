# ------------------------------------------------------------------------------
# Parse Olink proteomic data
# ------------------------------------------------------------------------------
# initiated on 2019-12-03
# authors :  Mun-Gwan Hong
# ------------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
# renv::install("Olink-Proteomics/OlinkRPackage/OlinkAnalyze")

#----- File names --------------------------------------------------------------
fn <- list(
  lib = "utils.R",
  i = list(                               #  input
    olink = list(
      INF = '../data/raw_internal/20210504/cederroth_tinnitus_protein_profiling_NPX_belowLOD.xlsx',
      NEU = '../data/raw_internal/20210505/Tinnitus2_Cederroth_NPX_below_LOD.xlsx'
    )
  ),
  o = list(                               #  output
    olk01 = "../data/s1-olink_proteomic.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(c(fn$lib, unlist(fn$i), dirname(unlist(fn$o))))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$olink$INF) == "42f834d68c672ff0bdbf5bfcc48fe68d",
  tools::md5sum(fn$i$olink$NEU) == "f8104973dfd7edf5d395d55b6a3a81c3"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

source(fn$lib) # as.matrix_olinkdf
#----- MAIN --------------------------------------------------------------------

# parse Olink data stored in two files
olink0 <- lapply(
  fn$i$olink, 
  OlinkAnalyze::read_NPX
) %>% 
  bind_rows()

i_ctrl <- unique(olink0$SampleID) %>% 
  str_subset("^Mix_|^IPC|^Neg Control")

olink1 <- olink0 %>% 
  mutate(MissingFreq = as.numeric(MissingFreq)) %>% # fix classes
  # extract plate ID and experiment date only
  separate(PlateID, into = c("PlateID", NA, NA, NA, "experiment_date"), sep = "_") %>% 
  mutate(
    PlateID = str_remove(PlateID, "^Plate") %>% as.integer,
    Panel = str_remove(Panel, "^Olink ") %>% 
      str_remove("Target 96 ") %>% 
      str_to_title()
  )

# Control samples
olink_ctrl <- olink1 %>% 
  filter(SampleID %in% i_ctrl)

# Olink without those controls
olink <- olink1 %>% 
  filter(!SampleID %in% i_ctrl)

# Save --------------------------------------------------------------------

save(olink, olink0, olink_ctrl, as.matrix_olinkdf, file = fn$o$olk01)
