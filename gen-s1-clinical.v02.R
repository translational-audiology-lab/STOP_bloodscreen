# -----------------------------------------------------------------------------#
# Limit the samples to those with proteomic data
# -----------------------------------------------------------------------------#
# initiated on 2021-05-21
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

#----- File names --------------------------------------------------------------

fn <- list(
  i = list(                               #  input
    olk02 = "../data/s1-olink_proteomic.v02.RData",
    c01 = "../data/s1-clinical.v01.RData"
  ),
  o = list(                               #  output
    c02 = "../data/s1-clinical.v02.RData"
  )
)
# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))

#----- MAIN --------------------------------------------------------------------

# Load proteomic data ver.02
load(fn$i$olk02)
# Load clinical information ver.01
load(fn$i$c01)

# confirm all samples in `olink` are also in `qns`
stopifnot(all(unique(olink$SampleID) %in% qns$RID))

qns <- qns %>% 
  # Limit to the 1088 samples with proteomic data
  filter(RID %in% unique(olink$SampleID)) %>% 
  # To avoid any confusion between two versions (1 and 2)
  rename(SampleID = RID)

# Save --------------------------------------------------------------------

save(qns, codes, about_table, file = fn$o$c02)
