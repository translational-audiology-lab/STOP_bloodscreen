# -----------------------------------------------------------------------------#
# Extract cleaned clinical info and proteomic data
# -----------------------------------------------------------------------------#
# initiated on 2021-05-30
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

#----- File names --------------------------------------------------------------

fn <- list(
  lib = "utils.R",
  i = list(                               #  input
    olk02 = "../data/s1-olink_proteomic.v02.RData",
    c02 = "../data/s1-clinical.v02.RData"
  ),
  o = list(                               #  output
    ut = "../data/exported_clean_data.xlsx"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(fn$lib, unlist(fn$i), dirname(unlist(fn$o)))))

#----- MAIN --------------------------------------------------------------------

# Load data
load(fn$i$olk02) # proteomic data 
load(fn$i$c02)  # clinical information
# confirm no missing samples
stopifnot(all(qns$SampleID %in% olink$SampleID))
stopifnot(all(olink$SampleID %in% qns$SampleID))

# Olink data in wide form
olink_wide <- olink %>% 
  unite(prot_olinkid, Assay, OlinkID, sep = " (") %>% 
  mutate(prot_olinkid = paste0(prot_olinkid, ")")) %>% 
  pivot_wider(SampleID, names_from = prot_olinkid, values_from = NPX)
# Olink wide form + clinical info
oc <- full_join(olink_wide, qns, by = "SampleID")

# Olink data only divided by panel
olink_by_panel <- olink %>% 
  split(.$Panel) %>% 
  lapply(
    . %>% 
      pivot_wider(SampleID, names_from = Assay, values_from = NPX)
  )

# save
out <- c(
  olink_by_panel,
  list("All data" = oc)
)
openxlsx::write.xlsx(out, file = fn$o$ut)
