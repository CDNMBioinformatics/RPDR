require(openxlsx)
require(data.table) # fread
require(dplyr) # left_join, mutate, mutate_at, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_plasma <- function(DF_to_fill = All_merged,
                           input_file_name = config$plasma_file_name,
                           date_cutoff,
                           reduce_to_plasma_only = TRUE){
  loginfo("Processing plasma dates file...")
  if (grepl("xlsx$", input_file_name)){
    plasma <- read.xlsx(input_file_name, detectDates = TRUE)
  } else {
    plasma <- fread(input_file_name)
  }
  if ("SUBJECTID" %in% names(plasma)){
    plasma <- plasma %>% rename("Biobank_Subject_ID" = SUBJECTID)
  }
  if ("coll_date" %in% names(plasma)){
    plasma <- plasma %>% rename("COLLECT_DATE" = coll_date)
  }
  plasma <- plasma %>% mutate(Biobank_Subject_ID = as.numeric(Biobank_Subject_ID),
                              COLLECT_DATE = ymd(COLLECT_DATE))
  if (!missing("date_cutoff")){
    plasma <- plasma %>% filter(COLLECT_DATE <= date_cutoff)
  }
  Output_Columns <- plasma %>% group_by(Biobank_Subject_ID) %>% arrange(COLLECT_DATE) %>%
    summarise(First_Collection_Date = first(COLLECT_DATE),
              Last_Collection_Date = last(COLLECT_DATE),
              nCollection_Dates = n(),
              All_Collection_Dates = paste(COLLECT_DATE, collapse = ";"),
              .groups = 'drop')
  if(reduce_to_plasma_only){
    DF_to_fill <- right_join(DF_to_fill, Output_Columns, by = "Biobank_Subject_ID")
  } else {
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "Biobank_Subject_ID")
  }
  return(DF_to_fill)
}
combine_and_process_plasma <- function(DF_to_fill = All_merged,
                                       before_file_name = config$before_plasma_file_name,
                                       after_file_name = config$after_plasma_file_name,
                                       date_cutoff,
                                       reduce_to_plasma_only = TRUE){
  loginfo("Processing plasma dates file...")
  if (grepl("xlsx$", before_file_name)){
    before_plasma <- read.xlsx(before_file_name, detectDates = TRUE)
  } else {
    before_plasma <- fread(before_file_name)
  }
  if (grepl("xlsx$", after_file_name)){
    after_plasma <- read.xlsx(after_file_name, detectDates = TRUE)
  } else {
    after_plasma <- fread(after_file_name)
  }
  if ("SUBJECTID" %in% names(before_plasma)){
    before_plasma <- before_plasma %>% rename("Biobank_Subject_ID" = SUBJECTID)
  }
  if ("SUBJECTID" %in% names(after_plasma)){
    after_plasma <- after_plasma %>% rename("Biobank_Subject_ID" = SUBJECTID)
  }
  if ("coll_date" %in% names(before_plasma)){
    before_plasma <- before_plasma %>% rename("COLLECT_DATE" = coll_date)
  }
  if ("coll_date" %in% names(after_plasma)){
    after_plasma <- after_plasma %>% rename("COLLECT_DATE" = coll_date)
  }
  plasma <- rbind(before_plasma, after_plasma)
  plasma <- unique(plasma)
  plasma <- plasma %>% mutate(Biobank_Subject_ID = as.numeric(Biobank_Subject_ID),
                              COLLECT_DATE = ymd(COLLECT_DATE))
  if (!missing("date_cutoff")){
    plasma <- plasma %>% filter(COLLECT_DATE <= date_cutoff)
  }
  Output_Columns <- plasma %>% group_by(Biobank_Subject_ID) %>% arrange(COLLECT_DATE) %>%
    summarise(First_Collection_Date = first(COLLECT_DATE),
              Last_Collection_Date = last(COLLECT_DATE),
              nCollection_Dates = n(),
              All_Collection_Dates = paste(COLLECT_DATE, collapse = ";"),
              .groups = 'drop')
  Output_Columns <- Output_Columns %>%
    mutate(Before_Plasma = ifelse(Biobank_Subject_ID %in% before_plasma$Biobank_Subject_ID,
                                  TRUE,
                                  FALSE),
           After_Plasma = ifelse(Biobank_Subject_ID %in% after_plasma$Biobank_Subject_ID,
                                 TRUE,
                                 FALSE),
           Any_Plasma = ifelse(Biobank_Subject_ID %in% before_plasma$Biobank_Subject_ID,
                               ifelse(Biobank_Subject_ID %in% after_plasma$Biobank_Subject_ID,
                                      "Both",
                                      "Before"),
                               "After"))
  if(reduce_to_plasma_only){
    DF_to_fill <- right_join(DF_to_fill, Output_Columns, by = "Biobank_Subject_ID")
  } else {
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "Biobank_Subject_ID")
  }
  return(DF_to_fill)
}