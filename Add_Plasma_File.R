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
  plasma <- read.xlsx(config$plasma_file_name, detectDates = TRUE)
  plasma <- plasma %>% rename("Biobank_Subject_ID" = SUBJECTID) %>%
    mutate(Biobank_Subject_ID = as.numeric(Biobank_Subject_ID),
           COLLECT_DATE = ymd(COLLECT_DATE))
  if (exists("date_cutoff")){
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