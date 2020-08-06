require(data.table) # fread
require(dplyr) # left_join, mutate, mutate_at, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

####################################################################################################
#################################### Deidentification functions ####################################
####################################################################################################
process_deidentified <- function(DF_to_fill = All_merged,
                                 input_file_name = config$biobank_file_name,
                                 PPV.NPV.only = FALSE,
                                 Asthma.only = FALSE,
                                 Asthma.COPD.only = FALSE,
                                 Asthma.Tobacco.only = FALSE,
                                 replace.existence.T.F = FALSE,
                                 clean.list = FALSE,
                                 get.date.range = FALSE){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name) %>%
    mutate(Biobank_Subject_ID = as.numeric(Biobank_Subject_ID)) # Just in case it reads as character
  # Remove non-word characters, remove duplicate, trailing, and leading spaces in names
  names(Deidentified) <- gsub("^_|_$", "", gsub("_+", "_", gsub("\\W", "_", names(Deidentified))))
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
  if (Asthma.COPD.only){
    loginfo("Select variables with Asthma or Tobacco use only")
    Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("Asthma"),
                                            contains("Tobacco"), contains("Smok"))
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
    return(List_df)
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
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "EMPI")
}