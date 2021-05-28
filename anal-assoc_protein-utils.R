# -----------------------------------------------------------------------------#
# Functions used in the association tests 
# -----------------------------------------------------------------------------#
# initiated on 2021-05-21
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

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
#' @param resample_t output of \code{\link{lm_prot_resample_t}}
#' @param olinkdf a data frame that contain Olink data in long format
#' @param c_data clinical data 
#' @param pos position of the variable of interest in the formula given in
#'   \code{rhs}. \code{pos = 1} is for the first independent variable.
#'
#' @return a tibble in which each row has output including the adjusted P-value
#'
#' @seealso [lm_out_1line()]
#' @references 
#' Westfall and Young's max-T method (1993)
#' \url{https://fdhidalgo.github.io/multitestr/articles/multitestr.html#stepdown}

lm_prot_padj_by_max_t <- function(resample_t, olinkdf, c_data, pos = 1) {
  keys <- attr(resample_t, "key")
  rhs <- keys$rhs
  
  # trim, because `c_data` can be huge
  c_data <- c_data %>% 
    select(SampleID, all_of(all.vars(formula(paste("~", rhs))))) %>% 
    droplevels()
  
  # make sure the same input was used in lm_prot_resample_t
  stopifnot(keys$c_data == digest::digest(c_data))
  
  # olink wide form + clinical info, common samples only
  oc <- olinkdf %>% 
    pivot_wider(SampleID, names_from = OlinkID, values_from = NPX) %>%   # wide
    inner_join(c_data, by = "SampleID")

  prots <- tibble(OlinkID = unique(olinkdf$OlinkID))    # proteins

  # observed T-statistics
  obs <- map_dfr(
    prots$OlinkID,
    function(ii) {
      formula(paste(ii, "~", rhs)) %>% 
        lm_out_1line(data = oc, pos = pos) %>% 
        mutate(OlinkID = ii, .before = 1L)
    }
  ) %>% 
    mutate(
      t_obs = abs(`t value`),
      bonf_P = p.adjust(Pval, method= 'bonferroni')
    )
  
  # max T during resampling
  Q_b <- resample_t[[pos]] %>%    # for one term only
    left_join(prots, ., by = "OlinkID") %>%   # limit to those in `olinkdf`
    column_to_rownames("OlinkID") %>%   # remove OlinkID column
    as.matrix() %>% 
    `[`(order(obs$t_obs), ) %>%    # from min T_obs to max T_obs
    abs() %>% 
    apply(2, cummax) %>% 
    `[`(rank(obs$t_obs), )     # return to original order
  
  obs %>% 
    mutate(
      # adjusted p-values by resampling
      P_perm = apply(Q_b >= t_obs, 1, mean) %>%
        `[`(order(t_obs, decreasing = TRUE)) %>%    # max T_obs -> min T_obs
        cummax() %>%          # successive maximization
        `[`(rank(-t_obs))         # return to original protein order
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
    select(..., "_p_ value" = Pval, all_of(after_p))
  
  # align numeric to right
  align <- sapply(res, is.numeric) %>% 
    if_else("r", "l") %>% 
    paste(collapse = "")
  
  res %>% 
    # # * 10^x format, but doesn't work with scroll
    # mutate(across(where(is.numeric), pval_toLatex)) %>% 
    mutate(across(where(is.numeric), 
                  ~ str_remove(format(.x, 3, digits = 3), "NA"))) %>% 
    kbl(caption = caption, align = align) %>% 
    kable_styling(full_width = F)
}








