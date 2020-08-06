require(data.table) # fread, fwrite
require(dplyr) # left_join, select, rename, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

## This file has the code that creates the first data frame of cleaned data

start_processing <- function(input_file_header = config$rpdr_file_header,
                             input_file_ending = config$rpdr_file_ending,
                             include.identifiable = FALSE){
  loginfo("Begin cleaning RPDR data...")
  loginfo("Processing demographics file...")
  Demographics <- fread(str_c(input_file_header, "Dem", input_file_ending))
  Demographics <- Demographics %>% select(EMPI, Gender, Race, Date_of_Birth, Age, Date_Of_Death, Vital_status)
  loginfo(str_c(nrow(Demographics), " subjects processed"))

  loginfo("Processing biobank ids file... ")
  BiobankIDs <- fread(str_c(input_file_header, "Bib", input_file_ending))
  BiobankIDs <- BiobankIDs %>% select(Subject_Id, EMPI) %>% rename(Biobank_Subject_ID = Subject_Id)
  loginfo(str_c(nrow(BiobankIDs), " subjects processed"))
  logwarning("if query was originally from RPDR, not all subjects may have Biobank Ids")
  
  loginfo("Generating merged output...")
  Merged_DF <- left_join(Demographics, BiobankIDs, by = "EMPI") %>%
    select(EMPI, Biobank_Subject_ID, everything())
  rm(Demographics, BiobankIDs)
  
  if(include.identifiable){
    loginfo("Processing identifiable information...")
    Identifiable <- fread(str_c(input_file_header, "Con", input_file_ending))
    Identifiable <- Identifiable %>% mutate(Zip = as.character(Zip),
                                            Zip = ifelse(grepl("CT|MA|ME|NH|NJ|RI|VT", State) &
                                                           (str_count(Zip) == 4 |  str_count(Zip) == 8),
                                                         str_c("0", Zip),
                                                         ifelse(grepl("PR|VI", State) & 
                                                                  (str_count(Zip) == 3 | str_count(Zip) == 7),
                                                                str_c("00", Zip),
                                                                Zip)),
                                            Zip = substr(Zip, 1, 5),
                                            Employee = ifelse(grepl("E", VIP), "YES", "NO"))
    Identifiable <- Identifiable %>% select(EMPI, Zip, Employee)
    loginfo(str_c(nrow(Identifiable), " subjects processed"))
    loginfo("Adding identified to merged output...")
    Merged_DF <- left_join(Merged_DF, Identifiable, by = "EMPI")
    rm(Identifiable)
  }
  loginfo("Merged output complete")
  return(Merged_DF)
}