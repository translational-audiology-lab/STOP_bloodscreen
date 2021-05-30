# -----------------------------------------------------------------------------#
# 1. Limit the samples to those with proteomic data
# 2. Add 'non_aids_analysis' for the analysis limited to non-hearing-aid users
# 3. `qset` questions sets (e.g. stress)
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


# With proteomic data -----------------------------------------------------


qns <- qns %>% 
  # Limit to the 1088 samples with proteomic data
  filter(RID %in% unique(olink$SampleID)) %>% 
  # To avoid any confusion between two versions (1 and 2)
  rename(SampleID = RID)


# Analysis of non-hearing-aid users ---------------------------------------
 
# Hearing aid users
aids_usr <- qns %>% 
  filter(A14_1 == "Yes") %>% 
  select(SampleID, Tinnitus, Sex, Age, Age_group, `Sample Lab`) %>% 
  # opposite tinnitus status for finding matched
  mutate(opp_tinnitus = fct_recode(Tinnitus, 
                                   "No" = "Yes, always", "Yes, always" = "No"))

print(dim(aids_usr))
print.data.frame(aids_usr %>% select(-opp_tinnitus))

# Matched pairs within hearing aid users
matched_pairs <- combn(1:nrow(aids_usr), 2) %>% 
  apply(2, function(ij) {
    if(aids_usr$Tinnitus[ij[1]] == aids_usr$opp_tinnitus[ij[2]] &
       aids_usr$Age[ij[1]] == aids_usr$Age[ij[2]] &
       aids_usr$Sex[ij[1]] == aids_usr$Sex[ij[2]]
    ) aids_usr$SampleID[c(ij[1], ij[2])] else NA
  }) %>% 
  subset(., !is.na(.))

# exclude pairs including duplicated IDs
tmp <- (which(duplicated(unlist(matched_pairs))) + 1) %/% 2
matched_pairs <- unlist(matched_pairs[-tmp])
print(matched_pairs)

# remaining hearing aid users
rem_aids_usr <- aids_usr %>% 
  filter(!SampleID %in% matched_pairs)

# Matched non-users to hearing users
set.seed(7)
matched_ids <- sapply(1:nrow(rem_aids_usr), function(ii) {
  qns %>% 
    select(A14_1, SampleID, Tinnitus, Sex, Age, Age_group, `Sample Lab`) %>% 
    filter(
      A14_1 == "No" &
        Tinnitus == rem_aids_usr$opp_tinnitus[ii] &
        Age == rem_aids_usr$Age[ii] &
        Sex == rem_aids_usr$Sex[ii]
    ) %>% 
    pull(SampleID) %>% 
    sample(1)
})
stopifnot(anyDuplicated(matched_ids) == 0)

# Mark for analysis limited to non-hearing aid users
qns <- qns %>% 
  mutate(
    non_aids_analysis = if_else(
      SampleID %in% c(aids_usr$SampleID, matched_ids),
      FALSE, TRUE
    )
  )


# Auestion sets -----------------------------------------------------------

qsets = list()

# stress related variables
qsets$stress <- tribble(
  ~lab,                 ~desc,
  "PSQ Total score",    "Stress",
  "HADS_A Total score", "Anxiety",
  "HADS_D Total score", "Depression",
  "HQ Total score",     "Hyperacusis",
  "A15_4",              "Temporomandibular joint pain",
  "A15_1",              "headache"
) %>% 
  mutate(with_q = paste0("`", lab, "`"))

# Find sub-type variables 
qsets$tinnitus_subtypes <- list(
  cat = c("TSCHQ_5", "TSCHQ_7", "TSCHQ_8", "TSCHQ_9", "TSCHQ_15", "TSCHQ_18", 
          "TSCHQ_22", "TSCHQ_24", "TSCHQ_27"),
  qty = c("TSCHQ_12", "TSCHQ_16", "TSCHQ_17")
)

# find Tinnitus Handicap Inventory scores
qsets$THI_scores <- Filter(function(x) grepl("^THI_", x), names(qns))


# Save --------------------------------------------------------------------

save(qns, codes, qsets, about_table, file = fn$o$c02)
