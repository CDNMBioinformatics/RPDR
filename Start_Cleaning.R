require(data.table) # fread, fwrite
require(dplyr) # left_join, select, rename, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

## This file has the code that creates the first data frame of cleaned data

start_processing <- function(input_file_header = config$rpdr_file_header,
                             input_file_ending = config$rpdr_file_ending,
                             include.identifiable = FALSE,
                             subset_ids){
  loginfo("Begin cleaning RPDR data...")
  loginfo("Processing demographics file...")
  Demographics <- fread(str_c(input_file_header, "Dem", input_file_ending))
  Demographics <- Demographics %>% select(EMPI, Gender, Race, Date_of_Birth, Age, Date_Of_Death, Vital_status) %>%
    mutate(Living = ifelse(grepl("Not reported as deceased", Vital_status), "Yes", "No"),
           Age_Range = case_when(Age < 10 ~ "0-9",
                                 Age < 20 ~ "10-19",
                                 Age < 30 ~ "20-29",
                                 Age < 40 ~ "30-39",
                                 Age < 50 ~ "40-49",
                                 Age < 60 ~ "50-59",
                                 Age < 70 ~ "60-69",
                                 Age < 80 ~ "70-79",
                                 Age < 90 ~ "80-89",
                                 Age >= 90 ~ "90+"))

  loginfo(str_c(nrow(Demographics), " subjects processed"))

  loginfo("Processing biobank ids file... ")
  BiobankIDs <- fread(str_c(input_file_header, "Bib", input_file_ending))
  BiobankIDs <- BiobankIDs %>% select(Subject_Id, EMPI) %>% rename(Biobank_Subject_ID = Subject_Id)
  loginfo(str_c(nrow(BiobankIDs), " subjects processed"))
  logwarn("if query was originally from RPDR, not all subjects may have Biobank Ids")
  
  loginfo("Generating merged output...")
  Merged_DF <- left_join(Demographics, BiobankIDs, by = "EMPI") %>%
    select(EMPI, Biobank_Subject_ID, everything())
  rm(Demographics, BiobankIDs)
  
  if(include.identifiable){
    loginfo("Processing identifiable information...")
    Identifiable <- fread(str_c(input_file_header, "Con", input_file_ending))
    Identifiable <- Identifiable %>%
      mutate(Zip = as.character(Zip),
             Zip = case_when(grepl("CT|MA|ME|NH|NJ|RI|VT", State) &
                               (str_count(Zip) == 4 |  str_count(Zip) == 8) ~ str("0", Zip),
                             grepl("PR|VI", State) & 
                               (str_count(Zip) == 3 | str_count(Zip) == 7) ~ str("00", Zip),
                             TRUE ~ Zip),
             Zip = substr(Zip, 1, 5),
             Employee = ifelse(grepl("E", VIP), "YES", "NO"))
    Identifiable <- Identifiable %>% select(EMPI, Zip, Employee)
    loginfo(str_c(nrow(Identifiable), " subjects processed"))
    loginfo("Adding identified to merged output...")
    Merged_DF <- left_join(Merged_DF, Identifiable, by = "EMPI")
    rm(Identifiable)
  }
  if(!missing(subset_ids)){
    loginfo(str_c("Selecting subset of ids from ", deparse(substitute(subset_ids))))
    Merged_DF <- Merged_DF %>% filter(Biobank_Subject_ID %in% subset_ids)
  } else {loginfo("Use full list of")}
  
  loginfo("Merged output complete")
  return(Merged_DF)
}