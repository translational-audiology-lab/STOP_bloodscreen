# -----------------------------------------------------------------------------#
# QC and preprocess - details in `report-QC-Olink_proteomic_data.Rmd`
# -----------------------------------------------------------------------------#
# initiated on 2021-05-12
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

library(tidyverse)

#----- File names --------------------------------------------------------------
fn <- list(
  i = list(                               #  input
    olk01 = "../data/s1-olink_proteomic.v01.RData", 
    c01 = "../data/s1-clinical.v01.RData"
  ),
  o = list(                               #  output
    olk02 = "../data/s1-olink_proteomic.v02.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(c(unlist(fn$i), dirname(unlist(fn$o))))))
#----- MAIN --------------------------------------------------------------------

# Load proteomic data ver.01
load(fn$i$olk01)
# To confirm sample IDs
load(fn$i$c01)

panels <- unique(olink$Panel)

# A table of QC warnings given to samples per panel
sample_warning <- olink %>% 
  distinct(SampleID, Panel, QC_Warning) %>% 
  pivot_wider(SampleID, names_from = Panel, values_from = QC_Warning)

# confirm no NA - Identical samples in both panels
stopifnot(all(!is.na(sample_warning)))

samples_in_olink_only <- setdiff(sample_warning$SampleID, qns$RID)

# * The samples without clinical information were excluded in the following analysis * #
olink <- olink %>% 
  filter(! SampleID %in% samples_in_olink_only)

# Confirm a LOD was given per Olink assay
olink %>%
  distinct(OlinkID, LOD, MissingFreq) %>% 
  {stopifnot(anyDuplicated(.$OlinkID) == 0)}

# Compute the below LLOD proportion, which is different from `MissingFreq`
# due to the samples with missing values
below_lod_prop <- olink %>% 
  drop_na(NPX) %>% {
    list(
      by_assay = group_by(., OlinkID),
      by_sample = group_by(., SampleID, Panel)
    )
  } %>% 
  lapply(. %>% summarise(below_lod_prop = sum(NPX <= LOD) / n(), .groups = "drop"))

# Assays having too many LLOD
i_assys_too_many_llod <- below_lod_prop$by_assay %>% 
  filter(below_lod_prop > 0.5) %>% 
  pull(OlinkID)

# the sample having too many LLOD
samples_too_many_llod <- below_lod_prop$by_sample %>% 
  filter(below_lod_prop > 0.5)

# * QC below LLOD * #
olink <- olink %>%
  drop_na(NPX) %>% 
  # Exclude those assays with too many values below LLOD
  filter(! OlinkID %in% i_assys_too_many_llod) %>% 
  # Exclude those samples with too many values below LLOD
  left_join(samples_too_many_llod, by = c("SampleID", "Panel")) %>% 
  filter(is.na(below_lod_prop)) %>% 
  select(-below_lod_prop)

save(olink, olink_ctrl, as.matrix_olinkdf, file = fn$o$olk02)
