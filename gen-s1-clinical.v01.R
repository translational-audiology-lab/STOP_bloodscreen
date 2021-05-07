# -----------------------------------------------------------------------------#
# Generate a RData file that contains clinical data
# Some processing was included as below.
#   - Remove two duplicated samples and two matching controls
#   - Codebook, incl. fixing the code errors
#   - Too many NAs or identical or too few element
#   - Special variables : Tinnitus, Sex, Smoking
#   - Fix values < LLOD
#   - add BMI
# -----------------------------------------------------------------------------#
# initiated on 2019-12-03
# authors : Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

library(tidyverse)

#----- File names --------------------------------------------------------------
fn <- list(
  lib = "utils.R",
  i = list(                               #  input
    qns = '../data/raw_internal/20210504/All STOP questionnaire data 180118_v14_BloodAnalysis_v2.xlsx',
    # codebooks for the variables of the main table
    code1 = '../data/raw_internal/20210505/Variable Key STOP Frågeformulär LG 2018.xls',
    code2 = '../data/raw_internal/20210505/Variablekey_ESITSQ.xlsx'
  ),
  o = list(                               #  output
    c01 = "../data/s1-clinical.v01.RData"
  )
)

# Brief check if all files exist
stopifnot(all(file.exists(c(unlist(fn$i), fn$lib, dirname(unlist(fn$o))))))
# Warn any difference in input file(s) using MD5 hash
if(!all(
  tools::md5sum(fn$i$qns) == "0a044eb9027c2a8b9a55253d0c554f0e",
  tools::md5sum(fn$i$code1) == "ac46ce6ad999e84cb4d6749643c9ffe8",
  tools::md5sum(fn$i$code2) == "84ac66b6f8d76c3b4e441052ee115f5f"
)) warning(
  "One input file or more is not exactly identical to the one used during ",
  "development of this script."
)

source(fn$lib) # irregular_cols
#----- MAIN --------------------------------------------------------------------

# Codebooks ---------------------------------------------------------------

# Read code books
code1 <- readxl::read_xls(fn$i$code1)
code2 <- readxl::read_xlsx(fn$i$code2)

# Clean up
codes <- bind_rows(code1, code2) %>% 
  rename(Variable = `Variable Name`, Codes = `Value Codes`) %>% 
  mutate(
    # Fix special characters
    Codes_tbl = Codes %>% 
      # '\n' between characters e.g. 'Very\ndissatisfie' -> 'Very dissatisfie'
      gsub("([[:alpha:]])\n([[:alpha:]])", "\\1 \\2", .) %>% 
      str_split("\n") %>%  # separate individual codes
      lapply(
        function(x) {
          if(isTRUE(all.equal(x, "none"))) return(NA)
          x %>% 
            trimws(whitespace = "[ \t\r\n\u00a0]") %>%    # \u00a0 = non breakable space
            tibble(Code = .) %>% 
            separate(Code, c("Code", "Level"), sep= ' = ')
        }
      ),
    Label = Label %>% 
      #  remove long repeated texts in questionnaires for a group of questions
      sub('FTQ - Fear of Tinnitus QuestionnaireThis questionnaire will help us understand how you think and feel about your tinnitus condition. It enables us to examine how tinnitus affects you, what effect is has on your mood, your behavior, your attitude. Below you will find 17 statements. Please check the box next to each statement that you think applies to your current situation. - ', '', .) %>% 
      sub('Hyperacusis QuestionnaireIn the following questionnaire, check the box corresponding to the answer which best applies to you. - ', '', .) %>% 
      sub('Tinnitus Handicap InventoryThe purpose of this questionnaire is to identify difficulties that you may be experiencing because of your tinnitus. Please answer every question. Please do not skip any question. - ', '', .) %>% 
      
      #  strange leading '  - '
      trimws(whitespace = "[ \t\r\n\u00a0]") %>% 
      sub("^-[ \t\r\n\u00a0]", '', .)
  )

# Main data table ---------------------------------------------------------

# Read the main data table including answers for questionnaires
qns0 <- readxl::read_xlsx(fn$i$qns, guess_max= 2000)
print(dim(qns0))

# Collect info about the given table through clean-up of the table
about_table <- list()

qns <- qns0 %>% 
  # Remove two duplicated samples and two matching controls 
  # (email : 20191125 from Christopher)
  filter(
    ! ID %in% c('US680082', 'VB209204',
                'EF746947', 'AJ899088')
  ) %>% 

  # remove NPX values.
  select(-c(IL8:`CSF-1`, `Plate ID`, `QC Warning`)) %>%
  
  mutate( # ID as character, Date as date, Time as time
    across(c('Barcode', 'ParticipantID'), as.character),   # IDs
    `Collection Date` = strptime(`Collection Date`, format= "%Y-%m-%dT%H:%M") %>% 
      as.POSIXct(),
    `Collection Time` = hms::as_hms(`Collection Date`)
  )

## Remove redundant ID variables ----------

# find them
about_table$redun_ids <- c(
  "ID of All STOP questionnaire data 180118_v6",
  "Stop_ID",
  "ID of All STOP questionnaire data 180118_v11",
  "ID of Untitled",
  "Barcode",
  "Assay"          
)

# Check a few redundant columns that have basically the same data
stopifnot(
  # Confirm `Assay` and `RID` were identical except three NAs in Assay
  isTRUE(
    all.equal(qns$RID[!is.na(qns$Assay)], qns$Assay[!is.na(qns$Assay)])
  ),
  identical(qns[["ID"]], qns[["Stop_ID"]])
)

cat("ID variables to be removed \n")
print(about_table$redun_ids)

qns <- qns %>% 
  select(-one_of(about_table$redun_ids))  # remove the redundant IDs

## Irregular columns ----------

about_table <- c(about_table, irregular_cols(qns, na_cuff_off = 0.8))

cat("> about_table\n")
print(names(about_table))

# Remove the columns having 
#  1) Mostly NA (>80%), including NA only
#  2) one identical value for all entries
#  3) one identical value or NA for all entries
qns <- qns %>% 
  select(-one_of(about_table$mostly_na, 
                 about_table$one_value_only$Column,
                 about_table$one_value_or_na$Column)) 

## Update with the codebook ----------

# Add the description of the codes from the codebook to the main table by
# changing variables from code to the meaning of the code
# - while fixing wrong codes in the codebook

#  Write down about wrong codes
about_table$wrong_code <- list(
  #  Yes or No wrong code (2, 1) instead of (1, 0)
  noyes = character(),
  #  'No', 'Yes, a little', 'Yes, quite a lot', 'Yes, a lot'
  #  (1, 2, 3, 4) instead of used (0, 1, 2, 3)
  noyes_little_to_lot = character()
)

# store current status
qns_tmp <- qns 

for(ii in c(1:nrow(codes))[!is.na(codes$Codes_tbl)]) { # for each variable in codebook
  if(codes$Variable[ii] %in% names(qns)) {
    this <- as.character(qns[[codes$Variable[ii]]])
    
    if(! all(this %in% codes$Codes_tbl[[ii]]$Code)) {
      #  Handle the code for missing value
      this[this == codes$`Missing Code`[ii]] <- NA_character_
      
      ##  Handle the errors in the codebook
      if(! all(this[!is.na(this)] %in% codes$Codes_tbl[[ii]]$Code)) {
        if(identical(codes$Codes_tbl[[ii]]$Level, c('Yes', 'No'))) {
          about_table$wrong_code$noyes <- c(about_table$wrong_code$noyes, codes$Variable[ii])
          codes$Codes_tbl[[ii]]$Code <- c('1', '0')
        }
        if(identical(codes$Codes_tbl[[ii]]$Level, 
                     c('No', 'Yes, a little', 'Yes, quite a lot', 'Yes, a lot'))) {
          about_table$wrong_code$noyes_little_to_lot <- c(
            about_table$wrong_code$noyes_little_to_lot, codes$Variable[ii]
          )
          codes$Codes_tbl[[ii]]$Code <- c('0', '1', '2', '3')
        }
        
        #  Confirm no weird code that was in the list of the codebook.
        stopifnot(all(this[!is.na(this)] %in% codes$Codes_tbl[[ii]]$Code))
      }
    }
    
    #  Change to factor
    this <- factor(
      this,
      levels = codes$Codes_tbl[[ii]]$Code,
      labels = codes$Codes_tbl[[ii]]$Level
    )
    
    #  Mark "Don't know' as missing
    levels(this)[tolower(levels(this)) == "don't know"] <- NA
    levels(this)[tolower(levels(this)) == "do not know"] <- NA
    
    #  Change to ordered factor
    if(
      codes$Variable[ii] %in% c("Intro_8", "A12", "A13", "A17", "B1", "B4") ||
      
      identical(levels(this), c('No', 'Yes, a little', 'Yes, quite a lot', 'Yes, a lot')) ||
      identical(levels(this), c("Less than 1 year", 
                                "1 to 2 years",
                                "2 to 5 years", 
                                "5 to 10 years", 
                                "More than 10 years")) ||
      identical(levels(this)[c(1, 3, 5)], c("Never", "Sometimes", "Always")) ||
      identical(levels(this)[1:3], c("Almost Never", "Sometimes", "Often")) ||
      #  THI-??
      identical(levels(this)[1:3], c("Yes", "Sometimes", "No")) ||
      
      identical(levels(this)[1:3], c("Not at all", "A little", "A moderate amount")) ||
      identical(levels(this)[1:3], c("Not at all", "A little", "Moderately")) ||
      identical(levels(this)[c(1, 2, 4)], c("Very dissatisfied", 
                                            "Dissatisfied",
                                            "Satisfied")) ||
      identical(levels(this)[2:10], format(1:9))
      
    ) {
      this <- ordered(this, levels= levels(this))
    } else {
      this <- droplevels(this)
    }
    
    qns[[codes$Variable[ii]]] <- this
  }
}

## Rename and fix Tinnitus, Sex, Smoking ----------

# confirm no `Tinnitus`, 'Sex', 'Smoking' before adding the variables
stopifnot(all(!c('Tinnitus', "Sex", "Smoking") %in% names(qns_tmp)))

qns <- qns %>% 
  rename(Smoking = A7, Sex = TSCHQ_2) %>% 
  
  mutate(
    #  Tinnitus status
    Tinnitus = Intro_3 %>% recode("Yes, always (all the time)" = "Yes, always")
  ) %>% 
  select(-Intro_3)

codes <- codes %>% 
  mutate(
    Variable = Variable %>% 
      replace(. == "A7", 'Smoking') %>% 
      replace(. == "TSCHQ_2", 'Sex') %>% 
      replace(. == "Intro_3", 'Tinnitus')
  )

## Derive additional variable(s) ----------

# BMI
qns <- qns %>% 
  mutate(
    BMI = A4 / (A3 / 100)^2
  )

## Clean-up more ----------
qns <- qns %>% 
  mutate(
    # Merge `MEB plan 4` with `Stockholm`
    `Sample Lab` = if_else(`Sample Lab` == 'MEB Plan 4', 'Stockholm', `Sample Lab`)
  )

## Remove the variables about existence of samples
has_vars <- names(qns) %>% {.[startsWith(., "Has")]}
print(has_vars)
qns <- qns %>% 
  select(-one_of(has_vars))

## Date/Time variables
datetime_vars <- names(qns) %>% 
  {.[grepl("time|date", tolower(.))]} %>% 
  {.[. != "Date difference"]}
qns <- qns %>% 
  mutate(
    across(
      all_of(datetime_vars),
      function(.x) {
        if(lubridate::is.POSIXt(.x) | lubridate::is.Date(.x) | hms::is_hms(.x)) return(.x)
        if(nchar(.x[1]) > 10) as.POSIXct(.x) else as.Date(.x)
      }
    )
  )

## Fix nations
variations <- list(
  sverige = c(
    'Sverige', 
    'sverige', 
    'Sverige, Västernorrland', 
    'Stockholm', 
    'Sverv', 
    'SE', 
    'SVERIGE', 
    'Svergie', 
    'Sveruge', 
    'Sv', 
    'Sverioge', 
    'Sverie', 
    'Sverige ä', 
    'Sveroge', 
    'Sverige, Stockholm', 
    'Svweige', 
    'Sverige och ca 3 mån per år i USA', 
    'Sveige', 
    'Scerige', 
    'Sverige/USA', 
    'Sverige 6 månader och Algarve, Portugal 6 månader', 
    "Kosovo sedan sex manader   (Normalt i Sverige)",
    "Göteborg",
    "Nyköping",
    "Kiruna"
  )
)
qns <- qns %>% 
  mutate_at(
    c("Intro_2", "Intro_2B"),
    function(.x) {
      .x[.x %in% variations$sverige] <- variations$sverige[1]
      .x[tolower(.x) == "usa"] <- "USA"
      .x[tolower(.x) == "iran"] <- "Iran"
      .x[.x == "Filand"] <- "Finland"
      .x[.x == "Storbrittanien"] <- "Storbritannien"
      .x[.x == "Prag"] <- "Tjeckoslovakien"
      .x[.x == "bosnien hercegovina"] <- "Bosnien"
      .x
    }
  )

# Save --------------------------------------------------------------------

save(qns, codes, about_table, file = fn$o$c01)


