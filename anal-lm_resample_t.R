# -----------------------------------------------------------------------------#
# Collect T statistics from linear regression models by resampling
# -----------------------------------------------------------------------------#
# initiated on 2021-05-25
# authors : Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

## CONSTANT
N_PERMUTATION = 10000


## Pre-loaded packages
# Please note that some functions are called directly from the package using
# double-colon '::' or triple-colon ':::' to make the source clearer.
library(tidyverse)
library(foreach)  # %dopar% foreach
library(doParallel)

## File names
fn <- list(
  i = list(                               # input
    olk02 = "../data/s1-olink_proteomic.v02.RData",
    c02 = "../data/s1-clinical.v02.RData"
  ),
  o = list(                               #  output
    perm_t = list(
      all_panels = "../cache/lm_resample_t.RData",
      Inflammation = "../cache/lm_resample_t-Inflammation.RData",
      Neurology = "../cache/lm_resample_t-Neurology.RData"
    )
  )
)
# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))


# Functions for resampling ------------------------------------------------

# ------------------
#' Get T stat only from [lm()]
#' 
#' This is faster than `coef(summary(lm()))` by skipping unnecessary steps and
#' output
#'
#' @param formula,data same as [lm()]
#' @return a vector of T values
#' 
#' @examples
#' \dontrun{
#' microbenchmark::microbenchmark(
#'   lm = coef(summary(lm(Sepal.Length ~ Species, iris)))[, "t value"],
#'   lm_t_values = lm_t_values(Sepal.Length ~ Species, iris),
#'   check = 'equal'
#' )
#' }
# ------------------

lm_t_values <- function(formula, data) {
  # decompose `lm` to speed up
  mf <- model.frame(formula, data, na.action = na.omit)
  y <- model.response(mf, "numeric")
  x <- model.matrix(formula, mf)
  z <- lm.fit(x = x, y = y)
  
  # decompose `summary.lm` to speed up
  p <- z$rank
  r <- z$residuals
  rss <- sum(r^2)
  rdf <- z$df.residual
  resvar <- rss/rdf
  Qr <- qr(x)
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients
  est/se
}


# ------------------
#' Collect T from resampling
#' 
#' This is for adjusting P-value by max-T resampling-based step-down procedure.
#' Because [purrr::nest()]is much slower than wide-form, [tidyr::pivot_wider()]
#' was chosen.
#'
#' @param rhs Right hand side of the formula. Left hand side is fixed to "NPX".
#' @param olinkdf a data frame that contain Olink data in long format
#' @param c_data clinical data 
#' @param count number of iterations
#'
#' @return the tibble of T statistics for the first term in `rhs` from
#'   resampling. The first column of the tibble has the OlinkIDs. This is to be
#'   used by [lm_prot_padj_by_max_t()]
#'
#' @import parallel foreach doParallel iterators digest
#'
#' @seealso [lm_prot_padj_by_max_t()]
#' @references 
#' Westfall and Young's max-T method (1993)
#' \url{https://fdhidalgo.github.io/multitestr/articles/multitestr.html#stepdown}
# ------------------

lm_prot_resample_t <- function(rhs, olinkdf, c_data, count = NULL) {
  # make sure unique sample IDs in the clinical data table
  stopifnot(anyDuplicated(c_data$SampleID) == 0)
  cat("~", rhs, "\n")  # progress

  # wide form
  olinkdf_wide0 <- olinkdf %>% 
    drop_na(NPX) %>% 
    pivot_wider(SampleID, names_from = OlinkID, values_from = NPX) 
  
  # trim, because `c_data` can be hugely wide
  terms_in_rhs <- all.vars(formula(paste("~", rhs)))
  c_data <- c_data %>% 
    # at least one protein data exists
    filter(SampleID %in% olinkdf_wide0$SampleID) %>% 
    select(SampleID, all_of(terms_in_rhs)) %>% 
    drop_na() %>%      # only complete cases
    droplevels()

  # clinical info complete cases only
  olinkdf_wide <- olinkdf_wide0 %>% 
    filter(SampleID %in% c_data$SampleID) %>% 
    select(-SampleID)   # no need of SampleID
  
  # proteins
  prots <- setNames(nm = unique(olinkdf$OlinkID))

  # backend for parallel computing
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)

  # T-statistics collected during resampling
  T_rnd <- foreach(
    # resampling
    rnd = iterators::isample(1:nrow(c_data), count = count),
    .combine = 'cbind',
    .inorder = FALSE,
    .packages = c("dplyr", "purrr")
  ) %dopar% {
    # Olink + randomized clinical data
    oc <- bind_cols(c_data[rnd, ], olinkdf_wide)
    # T-statistics from linear regression, the 1st term only
    map_dbl(prots, ~ lm_t_values(as.formula(paste(.x, "~", rhs)), oc)[2L])
  }

  # to make sure the same input
  attr(T_rnd, "key") <- list(
    rhs = rhs,
    c_data = digest::digest(c_data),
    olinkdf = digest::digest(olinkdf),
    dim = list(c_data = dim(c_data), olinkdf = dim(olinkdf))
  )
  T_rnd
}


# Load data ---------------------------------------------------------------

load(fn$i$olk02) # proteomic data 
load(fn$i$c02)  # clinical information
# confirm no samples in `qns` without proteomic data
stopifnot(all(qns$SampleID %in% olink$SampleID))


# Models to test ----------------------------------------------------------

stress_all <- paste(qsets$stress$with_q, collapse = " + ")

rhs <- c(
  "Tinnitus",
  "Tinnitus + Age + Sex + `Sample Lab`",
  "Tinnitus + Age + Sex + BMI",
  "Tinnitus + Age + Sex + Smoking",
  "Tinnitus + Age + Sex + BMI + `Sample Lab` + Smoking",
  "Tinnitus + Age + Sex + BMI + TSCHQ_26",
  "Tinnitus + Age + Sex + BMI + Smoking + TSCHQ_26",
  "Tinnitus + Age + Sex + BMI + `Sample Lab` + Smoking + TSCHQ_26",
  # Adjusted for stress related variables
  paste("Tinnitus + Age + Sex + BMI + Smoking + `Sample Lab` +", qsets$stress$with_q),
  paste("Tinnitus + Age + Sex + BMI + Smoking + `Sample Lab` +", stress_all),
  # hearing problem
  "TSCHQ_26 + Age + Sex + `Sample Lab`",
  "TSCHQ_26 + Age + Sex + BMI + Smoking + `Sample Lab`",
  # stress
  paste(qsets$stress$with_q, "+ Tinnitus + Age + Sex + BMI + Smoking + `Sample Lab`")
) %>% 
  setNames(nm = .)   # to keep the names after `map` or `lapply`


# Run by panels -----------------------------------------------------------

for(ipanel in names(fn$o$perm_t)) {
  if(ipanel == "all_panels") {
    olinkdf <- olink
    c_data <- qns
  } else {
    olinkdf <- olink %>% 
      filter(Panel == ipanel) 
    c_data <- qns %>% 
      filter(SampleID %in% olinkdf$SampleID)
  }
  # Initialize the storage of T stats from resampling
  resample_t_res <- list()

  # Overall samples
  cat("\n>>>", ipanel, ": Overall samples\n\n")
  resample_t_res$overall <- map(
    rhs,
    ~ lm_prot_resample_t(.x, olinkdf, c_data, count = N_PERMUTATION)
  )
  # Save
  save(resample_t_res, file = fn$o$perm_t[[ipanel]])
  
    
  # Non hearing-aid users
  cat("\n>>>", ipanel, ": Non hearing-aid users\n\n")
  qns_sub <- filter(c_data, non_aids_analysis)
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  resample_t_res$non_aid_user <- map(
    setNames(nm = c(
      "Tinnitus + Age + Sex + `Sample Lab`",
      "Tinnitus + Age + Sex + BMI + `Sample Lab` + Smoking + TSCHQ_26"
    )),
    ~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION)
  )
  
  # Female
  cat("\n>>>", ipanel, ": Female\n\n")
  qns_sub <- filter(c_data, Sex == "Female")
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  # remove " + Sex"
  resample_t_res$female <- setNames(nm = sub(" \\+ Sex", "", rhs)) %>% 
    map(~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION))
  
  # Male
  cat("\n>>>", ipanel, ": Male\n\n")
  qns_sub <- filter(c_data, Sex == "Male")
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  resample_t_res$male <- setNames(nm = sub(" \\+ Sex", "", rhs)) %>% 
    map(~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION))

  # Save
  save(resample_t_res, file = fn$o$perm_t[[ipanel]])

}

