require(data.table) # fread
require(dplyr) # left_join, mutate, mutate_at, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

#' Add deidentified columns from biobank file to dataset
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged).
#' @param input_file_name The path to the deidentified file (default = config$biobank_file_name).
#' @param PPV.NPV.only A logical value (default = FALSE), If TRUE only add columns with PPV or NPV in the name
#' @param Asthma.only A logical value (default = FALSE), If TRUE only add columns with Asthma in the name 
#' @param Asthma.COPD.only A logical value (default = FALSE), If TRUE only add columns with Asthma or COPD in the name
#' @param Asthma.Tobacco.only A logical value (default = FALSE), If TRUE only add columns with Asthma, Tobacco, Smoke, or Smoking in the name
#' @param Plasma.only A logical value (default = FALSE), If TRUE only add columns with Plasma in the name
#' @param replace.existence.T.F A logical value (default = FALSE), If TRUE replace Yes/No to 1/0
#' @param clean.list A logical value (default = FALSE), If TRUE update list format from [x1] [x2] to x1;x2
#' @param get.date.range A logical value (default = FALSE), If TRUE create year range columns based on pairs of columns with the same name but different dates
#' @param subset_ids NULL or a vector of biobank ids
#'
#' @return \code{DF_to_fill} modified with additional columns.
#' 
#' @examples
#' process_deidentified()
#' process_deidentified(Asthma.COPD.only = TRUE, clean.list = TRUE)
process_deidentified <- function(DF_to_fill = All_merged,
                                 input_file_name = config$biobank_file_name,
                                 PPV.NPV.only = FALSE,
                                 Asthma.only = FALSE,
                                 Asthma.COPD.only = FALSE,
                                 Asthma.Tobacco.only = FALSE,
                                 Plasma.only = FALSE,
                                 replace.existence.T.F = FALSE,
                                 clean.list = FALSE,
                                 get.date.range = FALSE,
                                 subset_ids){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  # Remove non-word characters, remove duplicate, trailing, and leading spaces in names
  names(Deidentified) <- gsub("^_|_$", "", gsub("_+", "_", gsub("\\W", "_", names(Deidentified))))
  # Just in case R reads column as characters instead of numerical
  Deidentified <- Deidentified %>% mutate(Biobank_Subject_ID = as.numeric(Biobank_Subject_ID))
  if(!missing(subset_ids)){
    loginfo(str_c("Selecting subset of ids from ", deparse(substitute(subset_ids))))
    Deidentified <- Deidentified %>% filter(Biobank_Subject_ID %in% subset_ids)
  }
  if (PPV.NPV.only){
    loginfo("Select variables with PPV/NPV only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("PPV"), contains("NPV"))
  }
  if (Asthma.only){
    loginfo("Select variables with Asthma only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("Asthma"))
  }
  if (Asthma.COPD.only){
    loginfo("Select variables with Asthma or COPD only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("Asthma"), contains("COPD"))
  }
  if (Asthma.Tobacco.only){
    loginfo("Select variables with Asthma or Tobacco use only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("Asthma"),
                                            contains("Tobacco"), contains("Smok"))
  }
  if (Plasma.only){
    loginfo("Select variables with plasma only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("Plasma"))
  }
  if (replace.existence.T.F){
    loginfo("Replace 'Yes/No' with 1/0")
    Deidentified <- Deidentified %>% 
      mutate_at(vars(contains("Existence")), ~ifelse(.x == "Yes", 1, 0))
  }
  if (clean.list){
    loginfo("Replaced bracked format 'List' columns to ; separated")
    Deidentified <- Deidentified %>%
      mutate_at(vars(contains("List")), ~gsub("\\[((\\d|\\.)+)\\]", "\\1", .x),
                vars(contains("List")), ~gsub(" ", ";", .x),
                vars(contains("List")), ~gsub(";$", "", .x),
                vars(contains("List")), ~gsub("^$", NA, .x))
  }
  if (get.date.range){
    loginfo("Create Range_Years for Date_First Date_Most_Recent pairs")
    Possible_Date_Columns <- grep("Date_(First|Most_Recent)", names(Deidentified), value = TRUE)
    Possible_Pairs <- gsub("(.*_)Date_(First|Most_Recent)", "\\1", Possible_Date_Columns)
    Date_Pairs <- Possible_Pairs[duplicated(Possible_Pairs)]
    for (Group in Date_Pairs){
      Deidentified <- Deidentified %>%
        mutate(!!(as.symbol(str_c(Group, "Range_Years"))) := 
                 time_length(interval(!!(as.symbol(str_c(Group, "Date_First"))),
                                      !!(as.symbol(str_c(Group, "Date_Most_Recent")))),
                             unit = "year")) %>% 
        select(Biobank_Subject_ID:matches(str_c("^", Group, "Date_First")),
               str_c(Group, "Date_Most_Recent"),
               str_c(Group, "Range_Years"), everything())
    }
    rm(Possible_Date_Columns, Possible_Pairs, Date_Pairs)
  }
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  return(DF_to_fill)
}