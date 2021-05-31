# -----------------------------------------------------------------------------#
# Functions used in reports including association tests 
# -----------------------------------------------------------------------------#
# initiated on 2021-05-21
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)


# Test results into one line ----------------------------------------------

# ------------------
#' Extract one line output
#' 
#' @describeIn out_1line Extract one line output from summary output
#'
#' @param formula,data formula and data frame as described for [lm()]
#' @param main_out Main output from the model in a numerical matrix
#' @param i_row extracting row, which contains results of the variable of
#'   interest
#' @param pval_c p-value column name. It will be replace with `Pval` to
#'   standardize it
#'
#' @return a tibble having one row, which has the results of the variable of
#'   interest only from the model

out_1line <- function(main_out, formula, i_row, pval_c) {
  main_out %>% 
    as_tibble(rownames = 'Term') %>% 
    slice(i_row) %>%
    mutate(
      Term = gsub("`", "", trimws(Term)),   # clean
      # formula to string
      Model = paste(as.character(formula)[c(2, 1, 3)], collapse = " ")
    ) %>%
    rename('Pval' = all_of(pval_c)) %>%     # standardize
    select(Model, Term, everything())
}

#' @describeIn out_1line Extract from linear regression
#'
#' @param pos the position of the variable of interest in the formula given in
#'   `formula`. The default `pos = 1L` is for the first independent variable.
#' @examples 
#' lm_out_1line(Sepal.Length ~ Species, iris)
#' lm_out_1line(Sepal.Length ~ Species, iris, 1:2)   # multiple lines are possible

lm_out_1line <- function(formula, data, pos = 1L) {
  coef(summary(lm(formula, data))) %>% 
    out_1line(formula, pos + 1, 'Pr(>|t|)')      # skip (intercept)
}

#' @describeIn out_1line Extract from ANOVA
#' @examples 
#' aov_out_1line(Sepal.Length ~ Species, iris)

aov_out_1line <- function(formula, data) {
  summary(aov(formula, data= data))[[1]] %>%  # output of `summary.aoc` is a list
    out_1line(formula, n() - 1, 'Pr(>F)')    # `aov` gives type-1 ANOVA results
}


#' @describeIn out_1line Extract from ANOVA and Bartlett P-value
#' @examples 
#' aov_bartlett_out_1line(Sepal.Length ~ Species, iris)

aov_bartlett_out_1line <- function(formula, data) {
  # formula for Bartlett test
  all_t <- all.vars(formula)
  f_b <- update(formula, paste("~", all_t[length(all_t)]))   # last term only
  aov_out_1line(formula, data) %>% 
    mutate(
      `Bartlett P` = tryCatch(    # at least 2 per group
        bartlett.test(f_b, data)$p.value,
        error = function(cond) NA_real_
      )
    )
}
  


# ------------------
#' Perform tests and collect 1line outputs
#'
#' @param olinkdf_plus a data frame that contain Olink data in long format plus
#'   other relevant variables
#' @param rhs the right hand side of a formula. The response variable was fixed 
#'   as `NPX`
#' @param fun_1line function that generate one line output, e.g.
#'   [lm_out_1line()]
#'
#' @return a tibble with `OlinkID`, `Protein`, `Panel`, columns from `fun_1line`
#'   including `Pval`, `qval`, and `P_bonf`
#' @import qvalue

olink_assoc_1lines <- function(olinkdf_plus, rhs, fun_1line) {
  f <- formula(paste("NPX ~", rhs))
  olinkdf_plus %>% 
    # trim out unnecessary columns
    select(OlinkID, Protein = Assay, Panel, all_of(all.vars(f))) %>%
    nest(data = -c(OlinkID, Protein, Panel)) %>% 
    add_column(map_dfr(.$data, ~ fun_1line(f, .x))) %>% 
    mutate(
      qval = tryCatch(
        qvalue::qvalue(Pval),
        # ERROR: maximum p-value is smaller than lambda range.
        error = function(cond) qvalue::qvalue_truncp(Pval)
      )$qvalues,
      P_bonf = p.adjust(Pval, method= 'bonferroni')
    ) %>% 
    select(-data) %>% 
    arrange(Pval)
}


# max T -------------------------------------------------------------------

# ------------------
#' Adjust P-value by max-T resampling-based step-down procedure
#' 
#' Because `purrr::nest`is much slower than `pivot_wider`, 'pivot_wider" was 
#' chosen. This was separate from [lm_prot_resample_t()] because resampling
#' takes too much time.
#'
#' @param resample_t output of [lm_prot_resample_t()]
#' @param olinkdf a data frame that contain Olink data in long format
#' @param c_data a data frame of clinical data 
#'
#' @return a tibble in which each row has output including the adjusted P-value
#'   by max-T
#'
#' @seealso [lm_out_1line()], [lm_prot_resample_t()]
#' @references 
#' Westfall and Young's max-T method (1993)
#' \url{https://fdhidalgo.github.io/multitestr/articles/multitestr.html#stepdown}

lm_prot_padj_by_max_t <- function(resample_t, olinkdf, c_data) {
  keys <- attr(resample_t, "key")
  rhs <- keys$rhs
  
  # trim, because `c_data` can be hugely wide
  terms_in_rhs <- all.vars(formula(paste("~", rhs)))
  c_data <- c_data %>% 
    # at least one protein data exists
    filter(SampleID %in% unique(olinkdf$SampleID)) %>% 
    select(SampleID, all_of(terms_in_rhs)) %>% 
    drop_na() %>%      # only complete cases
    droplevels()

  # make sure the same input was used in `lm_prot_resample_t`
  stopifnot(keys$c_data == digest::digest(c_data),  # after the trim
            keys$olinkdf == digest::digest(olinkdf))
  
  # proteins
  prots <- distinct(olinkdf, OlinkID, Protein = Assay, Panel)
  # Olink data in wide form
  olink_wide <- olinkdf %>% 
    pivot_wider(SampleID, names_from = OlinkID, values_from = NPX)
  # Olink wide form + clinical info, common samples only
  oc <- inner_join(olink_wide, c_data, by = "SampleID")

  # observed T-statistics
  obs <- map_dfr(
    colnames(olink_wide)[-1],    # proteins, excl. SampleID
    ~ formula(paste(.x, "~", rhs)) %>% 
      lm_out_1line(data = oc) %>%        # linear regression
      mutate(OlinkID = .x, .before = 1L)
  ) %>% 
    mutate(
      abs_t = abs(`t value`),
      P_bonf = p.adjust(Pval, method= 'bonferroni')
    ) %>% 
    right_join(prots, ., by = "OlinkID")  # add protein, assay info
  
  # max T during resampling
  stopifnot(identical(rownames(resample_t), obs$OlinkID))
  Q_b <- abs(resample_t) %>%
    `[`(order(obs$abs_t), ) %>%    # from min T_obs to max T_obs
    apply(2, cummax) %>% 
    `[`(rank(obs$abs_t), )     # return to original order
  
  # return
  obs %>% 
    mutate(
      # adjusted p-values by resampling
      P_perm = apply(Q_b >= abs_t, 1, mean) %>%   # more extreme proportion
        `[`(order(abs_t, decreasing = TRUE)) %>%    # max T_obs -> min T_obs
        cummax() %>%          # successive maximization
        `[`(rank(-abs_t))         # return to original protein order
    ) %>% 
    select(-abs_t)
}


# Show - out_df -----------------------------------------------------------

# ------------------
#' Show result and p values in a kable
#'
#' @param x a data frame that has `Pval`
#' @param ... a selection of column names in \code{x} to show
#' @param n number of hits to show. set `Inf` to show all
#' @param after_p columns to show after `p value`
#' @param caption caption in kable
#' 
#' @return A kable where rows were sorted by `Pval`

show_res_p <- function(
  x, ..., n = 3, caption = paste("Top", n, "proteins"), after_p = NULL
) {
  stopifnot("Pval" %in% names(x))
  
  res <- x %>%  
    slice_min(Pval, n = n) %>% 
    select(..., "_p_ value" = Pval, all_of(after_p))
  
  # align numeric to right
  align <- if_else(sapply(res, is.numeric), "r", "l") %>% 
    paste(collapse = "")
  
  res %>% 
    # # * 10^x format, but doesn't work with scroll
    # mutate(across(where(is.numeric), pval_toLatex)) %>% 
    mutate(across(where(is.numeric), 
                  ~ sub("NA", "", format(.x, 3, digits = 3)))) %>% 
    kbl(caption = caption, align = align) %>% 
    kable_styling(full_width = F)
}


# ------------------
#' Show the results of tinnitus association tests
#' 
#' Only for Tinnitus study
#'
#' @param out_df the data frame of the output from [lm_prot_padj_by_max_t()]
#' @param n,caption check [show_res_p()]
#' 
#' @seealso [lm_prot_padj_by_max_t()] [show_res_p()]
#' @import qvalue

show_delta_P_perm <- function(
  out_df, n = 5, caption = paste("Top", n, "proteins")
) {
  out_df %>% 
    mutate(
      qval = tryCatch(
        qvalue::qvalue(Pval),
        # ERROR: maximum p-value is smaller than lambda range.
        error = function(cond) qvalue::qvalue_truncp(Pval)
      )$qvalues
    ) %>% 
    show_res_p(
      Protein, 
      "$\\Delta$(Yes - No)" = Estimate, 
      "Adj.P by perm." = P_perm, 
      "_q_ value" = qval,
      n = n,
      caption = caption
    )
}


# ------------------
#' Show a summary table of the hits from regression models
#'
#' @param out_df a data frame to show. It is assumed to be generated by 
#'   [olink_assoc_1lines()]
#' @param ... arguments that will be passed to [kabelExtra::kbl()]
#'
#' @return the same as [kable()] output

show_hit_table <- function(out_df, ...) {
  x <- out_df %>% 
    select(-Model, -Term, -qval, -Panel) %>% 
    select(Protein, everything())
  
  names(x) <- names(x) %>% {
    case_when(
      # footnote for P_bonf
      . == "P_bonf" ~ paste0("P_bonf", footnote_marker_symbol(1)),
      TRUE ~ .
    )
  }

  dots <- list(...)
  if(is.null(dots$align)) {
    dots$align <- ifelse(sapply(x, is.numeric), "r", "l")
  }
  if(is.null(dots$escape)) dots$escape <- F   # footnote
  
  # # * 10^x format, but doesn't work with scroll
  # mutate(across(where(is.numeric), pval_toLatex)) %>% 
  dots$x <- mutate(x, across(where(is.numeric), ~ format(.x, 3, digits = 3)))
  do.call("kbl", dots) %>% 
    kable_styling(full_width = FALSE) %>% 
    footnote(symbol = "The `P_bonf` were Bonferroni adjusted P-values.")
}

# ------------------
#' Show a few examples in scatter plots
#'
#' @param oids Olink IDs
#' @param x a variable to show in X-axis
#' @param data given data frame

show_examples <- function(oids, x, data) {
  sorted_prots <- distinct(data, OlinkID, Assay) %>% 
    slice(match(oids, OlinkID)) %>% 
    pull(Assay)
  
  data %>% 
    filter(OlinkID %in% oids) %>% 
    drop_na({{x}}) %>% 
    mutate(Assay = factor(Assay, levels = sorted_prots)) %>%  # keep the order
    ggplot(aes(x = {{x}}, y = NPX)) +
    facet_wrap(~ Assay, scales = "free_y")
}


# RHS ---------------------------------------------------------------------

# ------------------
#' RHS to vector
#'
#' @param rhs the right hand side of a formula
#' @return a character vector, in which each has a term

rhs_to_vector <- function(rhs) {
  all.vars(formula(paste("~", rhs)))
}

# ------------------
#' Text for `fun` for protein
#'
#' @describeIn text_R_fun_prot Test in R
#' @param rhs the right hand side of a formula
#' @param fun function in text
#' @param ii the order of the highlighting term
#' @return character surrounded by backticks

text_R_fun_prot <- function(rhs, fun = "lm", ii = 1L) {
  stopifnot(length(ii) == 1L)
  rv <- rhs_to_vector(rhs)
  if(ii < 0) ii <- length(rv) + 1 + ii  # from the last
  # highlight with quotes
  rv[ii] <- paste0("'", rv[ii], "'")
  paste0("`", fun, "(Protein ~ ", paste(rv, collapse = " + "), ")`")
}

#' Only for Tinnitus study
#' 
#' @describeIn text_R_fun_prot Test in Rmd, exclude the first term
#' @param except the ith term
text_Rmd_terms <- function(rhs, except = 1) {
  rv <- rhs_to_vector(rhs)
  stopifnot(length(except) == 1L)
  if(except > 0) rv <- rv[-except]
  if(except < 0) rv <- rv[-(length(rv) + 1 + except)]  # from the last
  rv %>% 
    tolower() %>% {
      case_when(
        . == "tschq_26" ~ "hearing problem",
        . == "bmi" ~ "BMI",
        TRUE ~ .
      )
    } %>% 
    knitr::combine_words()
}






