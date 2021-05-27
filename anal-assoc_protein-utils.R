# -----------------------------------------------------------------------------#
# Functions used in the association tests 
# -----------------------------------------------------------------------------#
# initiated on 2021-05-21
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)
library(foreach)  # foreach, %dopar%

#' Extract one line output from summary output
#'
#' @param summary_out A matrix of `summary` output
#' @param f formula used for model
#' @param slice_r one row that contains results of the variable of interest
#' @param pval_c p-value column name
#'
#' @return a tibble of one row
out_1line <- function(summary_out, f, slice_r, pval_c) {
  summary_out %>% 
    as_tibble(rownames = 'Term') %>% 
    mutate(
      Term = trimws(Term) %>% 
        str_remove_all("`")
    ) %>% 
    slice(slice_r) %>%
    mutate("Model" = paste(f[2], f[3], sep= " ~ ")) %>%   # formula to string
    rename('Pval' = all_of(pval_c)) %>%    # standardize the name of p-value column
    select(Model, Term, everything())
}

#' Extract one line output from linear regression
#'
#' @param f,data formula and data frame that will be passed to \code{\link{lm}}
#' @param pos position of the variable of interest in the formula given in
#'   \code{f}. \code{pos = 1} is for the first independent variable.
#'
#' @return a tibble having one row, which is the row for the variable of interest
lm_out_1line <- function(f, data, pos = 1) {
  coef(summary(lm(f, data= data))) %>% 
    out_1line(f, pos + 1, 'Pr(>|t|)')      # skip (intercept)
}

#' Extract one line output from ANOVA
#'
#' @param f,data formula and data to be transferred to \code{\link{aov}}
#' @return a tibble having one row, which is the row for the variable of interest
aov_out_1line <- function(f, data) {
  summary(aov(f, data= data))[[1]] %>%  # output of `summary.aoc` is a list
    out_1line(f, n() - 1, 'Pr(>F)')    # `aov` gives type-1 ANOVA results
}


#' Extract T-statistic from `lm` summary
#'
#' @param f,data formula and data frame that will be passed to \code{\link{lm}}
#' @param pos position of the variable of interest in the formula given in
#'   \code{f}. \code{pos = 1} is for the first independent variable.
#'
#' @return T-statistic
lm_out_t_stat <- function(f, data, pos = 1) {
  lm(f, data= data) %>% 
    summary() %>% 
    coef() %>% 
    # as_tibble(rownames = 'Term') %>% 
    `[`(pos + 1, "t value")      # skip (intercept)
}


#' Adjust P-value by max-T resampling-based step-down procedure
#' 
#' Because `purrr::nest`is much slower than `pivot_wider`, 'pivot_wider" was 
#' chosen.
#'
#' @param rhs Right hand side of the formula that will be pass to
#'   \code{\link{lm_out_t_stat}} and \code{\link{lm_out_1line}}. Left hand side
#'   is fixed to "NPX".
#' @param olinkdf a data frame that contain Olink data in long format
#' @param c_data clinical data 
#' @param n_iter number of iterations
#' @param pos position of the variable of interest in the formula given in
#'   \code{rhs}. \code{pos = 1} is for the first independent variable.
#'
#' @return a tibble in which each row has output including the adjusted P-value
#' @import foreach
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators isample
#'
#' @seealso [lm_out_t_stat()], [lm_out_1line()]
#' @references 
#' Westfall and Young's max-T method (1993)
#' \url{https://fdhidalgo.github.io/multitestr/articles/multitestr.html#stepdown}

lm_prot_padj_by_max_t <- function(rhs, olinkdf, c_data, n_iter = 100, pos = 1) {
  # make sure unique sample IDs in the clinical data table
  stopifnot(anyDuplicated(c_data$SampleID) == 0)
  
  # proteins
  prots <- unique(olinkdf$OlinkID)
  # wide form
  olinkdf_wide <- olinkdf %>% 
    pivot_wider(SampleID, names_from = OlinkID, values_from = NPX)

  # observed T-statistics
  oc <- inner_join(olinkdf_wide, c_data, by = "SampleID")
  obs <- map_dfr(
    prots,
    function(ii) {
      f <- formula(paste(ii, "~", rhs))
      lm_out_1line(f, data = oc, pos = pos) %>% 
        mutate(OlinkID = ii, .before = Model)
    }
  ) %>% 
    mutate(
      t_obs = abs(`t value`),
      bonf_P = p.adjust(Pval, method= 'bonferroni')
    )

  # backend for parallel computing
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  
  # T-statistics collected during resampling
  T_rnd <- foreach(
    # resampling
    rnd = iterators::isample(c_data$SampleID, count = n_iter),
    .inorder = FALSE,
    .combine = 'cbind',
    .packages = c("dplyr", "purrr")
  ) %dopar% {
    # Olink + randomized clinical data
    oc <- c_data %>% 
      mutate(SampleID = rnd) %>% 
      inner_join(olinkdf_wide, ., by = "SampleID")
    
    # T-statistics from linear regression
    map_dbl(
      prots,
      function(ii) {
        f <- formula(paste(ii, "~", rhs))
        abs(lm_out_t_stat(f, oc, pos = pos))
      }
    )
  }

  # max T during resampling
  Q_b <- T_rnd[order(obs$t_obs), ] %>% 
    apply(2, cummax) %>% 
    `[`(rank(obs$t_obs), )
  
  # # `order` <-> `rank`
  # x <- c(sample(100), 20)
  # identical(x, x[order(x)][rank(x)])
  # identical(x, x[order(x, decreasing = T)][rank(-x)])
  
  obs %>% 
    mutate(
      # adjusted p-values by resampling
      perm_P = apply(Q_b >= t_obs, 1, function(x) sum(x) / length(x)) %>%
        # successive maximization
        {cummax(.[order(t_obs, decreasing = TRUE)])} %>% 
        # return to original protein order
        `[`(rank(-t_obs))
    ) %>% 
    select(-t_obs)
}


#' Show result and p values in a kable
#'
#' @param x a data frame that has `Pval`
#' @param ... a selection of column names in \code{x} to show
#' @param n number of hits to show. If \code{cutoff_q} is given, not used.
#' @param after_p columns to show after `p value`
#' @param caption caption in kable

show_res_p <- function(x, ..., n = 3, caption = paste("Top", n, "proteins"), after_p = NULL) {
  stopifnot("Pval" %in% names(x))

  res <- x %>%  
    slice_min(Pval, n = n) %>% 
    select(..., "p value" = Pval, all_of(after_p))
  
  # align numeric to right
  align <- sapply(res, is.numeric) %>% 
    if_else("r", "l") %>% 
    paste(collapse = "")
  
  res %>% 
    # # * 10^x format, but doesn't work with scroll
    # mutate(across(where(is.numeric), pval_toLatex)) %>% 
    mutate(across(where(is.numeric), ~ format(.x, 3, digits = 3))) %>% 
    kable(caption = caption, align = align) %>% 
    kable_styling(full_width = F)
}








