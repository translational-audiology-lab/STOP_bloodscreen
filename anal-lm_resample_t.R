# -----------------------------------------------------------------------------#
# Collect T statistics from linear regression models by resampling
# -----------------------------------------------------------------------------#
# initiated on 2021-05-25
# authors : Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

## Pre-loaded packages
# Please note that some functions are called directly from the package using
# double-colon '::' or triple-colon ':::' to make the source clearer.
library(tidyverse)
library(parallel)
library(foreach)  # %dopar%, foreach
library(doSNOW)
library(iterators)
library(progress)

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


# Function for resampling -------------------------------------------------

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
#' @return a list of the tibbles of T statistics in each resampling. Each
#'   element of the list is a term in the linear regression model. The first
#'   column of the tibble has the OlinkIDs.
#'
#' @import parallel foreach doSNOW iterators digest progress
#'
#' @references 
#' Westfall and Young's max-T method (1993)
#' \url{https://fdhidalgo.github.io/multitestr/articles/multitestr.html#stepdown}

lm_prot_resample_t <- function(rhs, olinkdf, c_data, count = NULL) {
  # make sure unique sample IDs in the clinical data table
  stopifnot(anyDuplicated(c_data$SampleID) == 0)
  
  # progress bar
  cat(">", rhs, "\n")
  pb <- txtProgressBar(max = count, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # wide form
  olinkdf_wide <- olinkdf %>% 
    pivot_wider(SampleID, names_from = OlinkID, values_from = NPX) %>% 
    select(-SampleID)   # no need of SampleID
  # proteins
  prots <- unique(olinkdf$OlinkID)
  
  # trim, because `c_data` can be huge
  c_data <- c_data %>% 
    select(SampleID, all_of(all.vars(formula(paste("~", rhs))))) %>% 
    droplevels()

  # backend for parallel computing
  n_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(n_cores)
  doSNOW::registerDoSNOW(cl)

  # T-statistics collected during resampling
  T_rnd <- foreach(
    # resampling
    rnd = iterators::isample(1:nrow(c_data), count = count),
    .inorder = FALSE,
    .packages = c("dplyr", "purrr"),
    .options.snow = opts   # progress bar
  ) %dopar% {
    # Olink + randomized clinical data
    oc <- bind_cols(c_data[rnd, ], olinkdf_wide)
    
    # T-statistics from linear regression
    map_dfr(
      prots,
      function(ii) {
        f <- formula(paste(ii, "~", rhs))
        
        # decompose `lm` to speed up
        mf <- model.frame(f, oc, na.action = na.omit)
        y <- model.response(mf, "numeric")
        x <- model.matrix(f, mf)
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
        (est/se)[-1] # skip (intercept)
      }
    )
  }
  # cluster close
  close(pb)
  parallel::stopCluster(cl)
  
  # column names to 'r0001', 'r0002', ...
  ndigit <- nchar(format(length(T_rnd)))
  sform <- paste0("r%0", ndigit, "d")
 
  terms_rhs <- names(T_rnd[[1]])   # terms in the model
  # Transform to list of terms, where each element is a tibble of T-stat. 
  # In the tibble, proteins are in rows and repeated runs are in columns.
  out <- lapply(seq(terms_rhs), function(ii) {
    map(T_rnd, ~ unname(.x[, ii])) %>% 
      bind_cols(.name_repair = "minimal") %>% # avoid new names message
      `colnames<-`(sprintf(sform, 1:length(T_rnd))) %>% 
      mutate(OlinkID = prots, .before = 1)
  }) %>% 
    `names<-`(terms_rhs)
  # to make sure the same input
  attr(out, "key") <- list(
    rhs = rhs,
    c_data = digest::digest(c_data)
  )
  out
}


# Load data ---------------------------------------------------------------

load(fn$i$olk02) # proteomic data 
load(fn$i$c02)  # clinical information
# confirm no samples in `qns` without proteomic data
stopifnot(all(qns$SampleID %in% olink$SampleID))

## CONSTANT
N_PERMUTATION = 100


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
  "TSCHQ_26 + Age + Sex + `Sample Lab`",
  "TSCHQ_26 + Age + Sex + BMI + Smoking + `Sample Lab`"
) %>% 
  `names<-`(., .)   # to keep the names after `map` or `lapply`


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
  
  # Non hearing-aid users
  cat("\n>>>", ipanel, ": Non hearing-aid users\n\n")
  qns_sub <- filter(c_data, non_aids_analysis)
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  resample_t_res$non_aid_user <- map(
    list(
      "Tinnitus + Age + Sex + `Sample Lab`",
      "Tinnitus + Age + Sex + BMI + `Sample Lab` + Smoking + TSCHQ_26"
    ) %>%
      `names<-`(., .),
    ~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION)
  )
  
  # Female
  cat("\n>>>", ipanel, ": Female\n\n")
  qns_sub <- filter(c_data, Sex == "Female")
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  resample_t_res$female <- rhs %>% 
    str_remove(" \\+ Sex") %>%   # remove " + Sex"
    `names<-`(., .) %>% 
    map(~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION))
  
  # Male
  cat("\n>>>", ipanel, ": Male\n\n")
  qns_sub <- filter(c_data, Sex == "Male")
  olink_sub <- filter(olinkdf, SampleID %in% qns_sub$SampleID)
  resample_t_res$male <- rhs %>% 
    str_remove(" \\+ Sex") %>%   # remove " + Sex"
    `names<-`(., .) %>% 
    map(~ lm_prot_resample_t(.x, olink_sub, qns_sub, count = N_PERMUTATION))

  # Save
  save(resample_t_res, file = fn$o$perm_t[[ipanel]])

}

