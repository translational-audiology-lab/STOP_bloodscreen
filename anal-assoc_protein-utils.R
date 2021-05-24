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

#' Adjust P-value by max-T resmapling-based step-down procedure
#'
#' @param f formula that will be pass to \code{\link{lm_out_t_stat}} and
#'   \code{\link{lm_out_1line}}
#' @param olinkdf a data frame that contain Olink data in long format
#' @param c_data clinical data 
#' @param n_iter number of iterations
#' @param pos position of the variable of interest in the formula given in
#'   \code{f}. \code{pos = 1} is for the first independent variable.
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

lm_padj_by_max_t <- function(f, olinkdf, c_data, n_iter = 100, pos = 1) {
  # make sure no change after 'left_join' with p_data
  stopifnot(anyDuplicated(c_data$SampleID) == 0)
  
  # minimize the data size
  olinkdf <- olinkdf %>% 
    select(SampleID, OlinkID, Assay, NPX)
  c_data <- c_data %>% 
    select(SampleID, all_of(all.vars(f)[all.vars(f) != "NPX"]))
  
  # observed T-statistics
  obs <- inner_join(olinkdf, c_data, by = "SampleID") %>%
    nest(data = -OlinkID) %>% 
    add_column(map_dfr(.$data, ~ lm_out_1line(f, .x))) %>%
    mutate(
      t_obs = abs(`t value`),
      bonf_P = p.adjust(Pval, method= 'bonferroni')
    ) %>% 
    select(-data)

  # backend for parallel computing
  doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
  
  # T-statistics collected during resampling
  T_rnd <- foreach(
    # resampling
    rnd = iterators::isample(c_data$SampleID, count = n_iter),
    .inorder = FALSE,
    .combine = 'cbind',
    .packages = c("dplyr", "tidyr", "purrr")
  ) %dopar% {
    # T-statistics from linear regression
    t_b <- c_data %>% 
      mutate(SampleID = rnd) %>% 
      left_join(olinkdf, ., by = "SampleID") %>% 
      nest(data = -c(OlinkID, Assay)) %>% 
      mutate(abs_t_rnd = abs(map_dbl(data, ~ lm_out_t_stat(f, .x))))
  
    stopifnot(identical(t_b$OlinkID, obs$OlinkID))
    t_b$abs_t_rnd
  }

  # max T during resampling
  Q_b = T_rnd[order(obs$t_obs), ] %>% 
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


#' Compute q values. If an error occurs, BH procedure is used.
#'
#' @param p a numeric vector of p values 
#'
#' @return an object of \code{\link{qvalue::qvalue}}
#' @importFrom qvalue qvalue
compute_qval <- function(p) {
  qvals <- try(qvalue::qvalue(p), silent= T)
  
  if(inherits(qvals, "try-error")) {
    qvals <- qvalue::qvalue(p, pi0 = 1)    # BH procedure
  }
  qvals$call <- match.call() # copy the call
  return(qvals)
}

#' Show p values and q values in a kable
#'
#' @param x a data frame that has `Pval`
#' @param ... a selection of column names in \code{x} to show
#' @param cutoff_q a cutoff based on q-value
#' @param n number of hits to show. If \code{cutoff_q} is given, not used.
#' @param caption caption in kable
#' @seealso [compute_pval()]
show_p_q <- function(x, ..., cutoff_q = NULL, n = 3, caption = paste("Top", n, "proteins")) {
  stopifnot("Pval" %in% names(x))
  if("qval" %in% names(x)) { # compute q-value unless exists already
    res <- x
  } else {
    qvals <- compute_qval(x$Pval)
    if(qvals$pi0 == 1) caption <- paste(caption, "(q value by BH)")
    qvals <- qvals$qvalues

    res <- x %>% 
      mutate(qval = qvals)
  }
  
  res <- if(is.null(cutoff_q)) {
    slice_min(res, Pval, n = n)
  } else {
    filter(res, qval < cutoff_q)
  }
  
  res %>% 
    arrange(Pval) %>% 
    select(..., "p value" = Pval, "q value" = qval) %>% 
    # # * 10^x format, but doesn't work with scroll
    # mutate(across(where(is.numeric), pval_toLatex)) %>% 
    mutate(across(where(is.numeric), ~ format(.x, 3, digits = 3))) %>% 
    kable(caption = caption) %>% 
    kable_styling(full_width = F)
}








