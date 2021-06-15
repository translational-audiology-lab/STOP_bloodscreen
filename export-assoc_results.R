# -----------------------------------------------------------------------------#
# Export association test results
# -----------------------------------------------------------------------------#
# initiated on 2021-06-15
# authors :  Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

# Please note some functions are called directly using '::' or ':::'
library(tidyverse)

#----- File names --------------------------------------------------------------
panels <- c("Inflammation", "Neurology") %>% 
  c(., paste(., collapse = "_"))
sex <- c("overall", "female", "male")
p_s <- paste(rep(panels, each = 3), sex, sep = "-")

fn <- list(
  i = list(                               #  input
    res = map(setNames(nm = p_s), ~ paste0("../cache/assoc-", .x, ".RData"))
  ),
  o = list(                               #  output
    ut = map(setNames(nm = p_s), 
             ~ paste0("../results/tables/assoc_test_results-", .x, ".xlsx"))
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(unlist(fn$i), dirname(unlist(fn$o)))))

#----- MAIN --------------------------------------------------------------------

#' Extract Term and Model from each output file
#' 
#' @description It searches the two result columns recursively. 
#'
#' @param x list or data frame of association test results
#' @return data frame having `Term` and `Model`
get_Model <- function(x) {
  if(is.data.frame(x)) {
    return(x %>% 
             select(Term, Model) %>% 
             slice(1) %>% 
             mutate(res = list(x))
    )
  }
  map_dfr(x, ~ get_Model(.x))
}


for(ii in p_s) {
  load(fn$i$res[[ii]])   # association test results
  
  # into a list
  sheets <- get_Model(assoc_res) %>% 
    mutate(
      Sheet = sprintf("S%02d", 1:n()),
      Model = str_replace(Model, "OID00471", "NPX"),
      # one table for all 25 THI_??
      Model = str_replace(Model, "THI_1", "THI_??"),
      Term = str_replace(Term, "THI_1.L", "THI_??.L")
    ) %>% 
    select(Sheet, Term, Model, res)
  
  # Save
  out <- c(
    list(Notes = sheets %>% select(-res)),
    setNames(sheets$res, sheets$Sheet)
  )
  openxlsx::write.xlsx(out, file = fn$o$ut[[ii]])
}
