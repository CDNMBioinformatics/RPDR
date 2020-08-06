require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_encounters <- function(DF_to_fill = All_merged,
                               input_file_header = config$rpdr_file_header,
                               input_file_ending = config$rpdr_file_ending,
                               path_enc_abn = str_c(config$data_dir, "Enc_abnormalities/"),
                               output_file_ending = config$general_file_ending){
  loginfo("Processing Encounters data...")
  Encounters <- data.table(fread(str_c(input_file_header, "Enc", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_enc_abn)) {dir.create(path_enc_abn)}
  logdebug(str_c("Note: All encounter abnormalites can be found at ", path_enc_abn))
  # Cleaning: reformat date objects and resort for easier visualization
  Encounters <- Encounters %>% mutate(Admit_Date = as.Date(Admit_Date, "%m/%d/%Y"),
                                      Discharge_Date = as.Date(Discharge_Date, "%m/%d/%Y")) %>%
    arrange(EMPI, desc(Discharge_Date), Admit_Date, desc(Encounter_Status))
  Output_columns <- Encounters %>% group_by(EMPI) %>% summarise(Encounter_Info_Available = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Encounter_Info_Available = ifelse(is.na(Encounter_Info_Available), "No", Encounter_Info_Available))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter and have a covid positive test"))
  
  Output_columns <- Encounters %>% filter(Inpatient_Outpatient == "Inpatient") %>% group_by(EMPI) %>%
    summarise(Inpatient_Encounter_Info_Available = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Inpatient_Encounter_Info_Available = ifelse(is.na(Inpatient_Encounter_Info_Available), "No", Inpatient_Encounter_Info_Available))
  loginfo(str_c(nrow(Output_columns), " subjects have had an inpatient encounter and have a covid positive test"))
  
  Output_columns <- Encounters %>% filter(Inpatient_Outpatient == "Outpatient") %>% group_by(EMPI) %>%
    summarise(Outpatient_Encounter_Info_Available = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Outpatient_Encounter_Info_Available = ifelse(is.na(Outpatient_Encounter_Info_Available), "No", Outpatient_Encounter_Info_Available))
  loginfo(str_c(nrow(Output_columns), " subjects have had an outpatient encounter and have a covid positive test"))
  
  Output_columns <- Encounters %>% filter(Inpatient_Outpatient == "Outpatient-Emergency") %>% group_by(EMPI) %>%
    summarise(Emergency_Encounter_Info_Available = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Emergency_Encounter_Info_Available = ifelse(is.na(Emergency_Encounter_Info_Available), "No", Emergency_Encounter_Info_Available))
  loginfo(str_c(nrow(Output_columns), " subjects have had an outpatient-emergency encounter and have a covid positive test"))
  
  # Cleaning: Let's assume that no covid diagnosis was made before January 01, 2020
  Encounters <- Encounters %>% filter(Admit_Date > "2020-01-01")
  Output_columns <- Encounters %>% group_by(EMPI) %>% summarise(Encounter_Info_After_01012020 = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Encounter_Info_After_01012020 = ifelse(is.na(Encounter_Info_After_01012020), "No", Encounter_Info_After_01012020))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter since Jan 01, 2020"))
  
  # Cleaning: We are looking for LOS of relevance so let's ignore encounter information
  # up to a week before the Encounter
  
  if (!grepl("Positive_Result", names(DF_to_fill))){
    source(str_c("config$shared_functions_dir", "Clean_Labs_File.R"))
    DF_to_fill <- process_Covid_labs(DF_to_fill = DF_to_fill)
  }
  
  Dates_of_Interest <- DF_to_fill %>% filter(Positive_Result == "Yes") %>%
    extract(First_Positive_Test_Date, c("First_Positive_Test_Date", "First_Positive_Test_Time"),
            regex = "(\\d{2}/\\d{2}/\\d{4}) (\\d{2}:\\d{2})", remove = FALSE) %>%
    select(EMPI, First_Positive_Test_Date) %>%
    mutate(First_Positive_Test_Date = as.Date(First_Positive_Test_Date, "%m/%d/%Y"),
           Lower_Bound_Date = First_Positive_Test_Date - 7)
  Encounters <- left_join(Encounters, Dates_of_Interest)
  Encounters <- Encounters %>% filter(Admit_Date >= Lower_Bound_Date)
  Output_columns <- Encounters %>% group_by(EMPI) %>% summarise(Encounter_Info_After_1_Week_Before_FirstTest = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Encounter_Info_After_1_Week_Before_FirstTest = ifelse(is.na(Encounter_Info_After_1_Week_Before_FirstTest), "No", Encounter_Info_After_1_Week_Before_FirstTest))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter after a week before their first positive test"))
  
  rm(Dates_of_Interest)
  
  # Cleaning: Let's only look at the lines that have something to do with COVID
  Encounters <- Encounters %>%
    mutate(contains_covid = grepl("COVID", Admitting_Diagnosis) | grepl("COVID", Principal_Diagnosis) |
             grepl("COVID", Diagnosis_1) | grepl("COVID", Diagnosis_2) |
             grepl("COVID", Diagnosis_3) | grepl("COVID", Diagnosis_4) |
             grepl("COVID", Diagnosis_5) | grepl("COVID", Diagnosis_6) |
             grepl("COVID", Diagnosis_7) | grepl("COVID", Diagnosis_8) |
             grepl("COVID", Diagnosis_9) | grepl("COVID", Diagnosis_10))
  Covid_Ids <- Encounters %>% filter(contains_covid) %>% group_by(EMPI) %>% summarise() %>% pull()
  Encounters <- Encounters %>% filter(EMPI %in% Covid_Ids)
  Output_columns <- Encounters %>% group_by(EMPI) %>% summarise(Covid_Related_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Related_Encounter_Info = ifelse(is.na(Covid_Related_Encounter_Info), "No", Covid_Related_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter that had a COVID diagnosis"))
  
  rm(Covid_Ids)
  
  # Compare Inpatient metrics
  Inpatient <- Encounters %>% filter(Inpatient_Outpatient == "Inpatient")
  Output_columns <- Inpatient %>% group_by(EMPI) %>% summarise(Covid_Related_with_Inpatient_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Related_with_Inpatient_Encounter_Info = ifelse(is.na(Covid_Related_with_Inpatient_Encounter_Info), "No", Covid_Related_with_Inpatient_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter that had a COVID diagnosis and at least one inpatient encounter"))
  
  Covid_Inpatient <- Inpatient %>% filter(contains_covid)
  Output_columns <- Covid_Inpatient %>% group_by(EMPI) %>% summarise(Covid_Specific_Inpatient_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Inpatient_Encounter_Info = ifelse(is.na(Covid_Specific_Inpatient_Encounter_Info), "No", Covid_Specific_Inpatient_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an inpatient encounter"))
  
  Covid_Inpatient_Status <- Covid_Inpatient %>% filter(LOS_Days > 0) %>% arrange(Admit_Date)
  Output_columns <- Covid_Inpatient_Status %>% group_by(EMPI) %>%
    summarise(Covid_Specific_Inpatient_Encounter_LOS_Information = "Yes",
              Covid_Specific_Inpatient_Encounter_LOS_Occurances = n(),
              Covid_Specific_Inpatient_Encounter_LOS_First = first(LOS_Days),
              Covid_Specific_Inpatient_Encounter_LOS_Last = last(LOS_Days),
              Covid_Specific_Inpatient_Encounter_LOS_Total = sum(LOS_Days),
              Covid_Specific_Inpatient_Encounter_LOS_First_Admit_to_Last_Discharge = last(Discharge_Date) - first(Admit_Date),
              Covid_Specific_Inpatient_Encounter_LOS_Time_in_between_First_and_Second_Stay = ifelse(Covid_Specific_Inpatient_Encounter_LOS_Occurances > 1,
                                                                                                    nth(Admit_Date, 2) - first(Discharge_Date),
                                                                                                    NA),
              Covid_Specific_Inpatient_Encounter_LOS_Time_in_between_Second_and_Third_Stay = ifelse(Covid_Specific_Inpatient_Encounter_LOS_Occurances > 2,
                                                                                                    nth(Admit_Date, 3) - nth(Discharge_Date, 2),
                                                                                                    NA),
              .groups = 'drop')
  if (sum(is.na(Output_columns$Covid_Specific_Inpatient_Encounter_LOS_Time_in_between_Second_and_Third_Stay))){
    Output_columns <- Output_columns %>% select(-Covid_Specific_Inpatient_Encounter_LOS_Time_in_between_Second_and_Third_Stay)
  }
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Inpatient_Encounter_LOS_Information = ifelse(is.na(Covid_Specific_Inpatient_Encounter_LOS_Information), "No", Covid_Specific_Inpatient_Encounter_LOS_Information))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an inpatient encounter and have LOS information"))
  
  rm(Inpatient, Covid_Inpatient, Covid_Inpatient_Status)
  
  # Compare Outpatient metrics
  Outpatient <- Encounters %>% filter(Inpatient_Outpatient == "Outpatient")
  Output_columns <- Outpatient %>% group_by(EMPI) %>% summarise(Covid_Related_with_Outpatient_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Related_with_Outpatient_Encounter_Info = ifelse(is.na(Covid_Related_with_Outpatient_Encounter_Info), "No", Covid_Related_with_Outpatient_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter that had a COVID diagnosis and at least one outpatient encounter"))
  
  Covid_Outpatient <- Outpatient %>% filter(contains_covid)
  Output_columns <- Covid_Outpatient %>% group_by(EMPI) %>% summarise(Covid_Specific_Outpatient_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Outpatient_Encounter_Info = ifelse(is.na(Covid_Specific_Outpatient_Encounter_Info), "No", Covid_Specific_Outpatient_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an outpatient encounter"))
  
  Covid_Outpatient_Status <- Covid_Outpatient %>% filter(Encounter_Status == "Regular" |
                                                           Encounter_Status == "Inhouse" | 
                                                           Encounter_Status == "Contains Invalid DRGs") %>%
    arrange(EMPI, Admit_Date, desc(Discharge_Date))
  Covid_Outpatient_Status <- distinct(Covid_Outpatient_Status, EMPI, Admit_Date, .keep_all = TRUE)
  Covid_Outpatient_Status <- distinct(Covid_Outpatient_Status, EMPI, Discharge_Date, .keep_all = TRUE)
  Output_columns <- Covid_Outpatient_Status %>% group_by(EMPI) %>%
    summarise(Covid_Specific_Outpatient_Encounter_Existence = "Yes",
              Covid_Specific_Outpatient_Encounter_Occurances = n(),
              Covid_Specific_Outpatient_Consecutive_Encounters_Existence = ifelse(any(ifelse(any(LOS_Days > 0), TRUE, FALSE),
                                                                                      ifelse(Covid_Specific_Outpatient_Encounter_Occurances > 1,
                                                                                             any(Admit_Date[2:length(Admit_Date)] - Discharge_Date[1:length(Discharge_Date) - 1] == 1),
                                                                                             FALSE)),
                                                                                  "Yes", "No"),
              Covid_Specific_Outpatient_LOS_Total = sum(LOS_Days) + ifelse(Covid_Specific_Outpatient_Encounter_Occurances > 1,
                                                                           sum(Admit_Date[2:length(Admit_Date)] - Discharge_Date[1:length(Discharge_Date) - 1] == 1),
                                                                           0),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Outpatient_Encounter_Existence = ifelse(is.na(Covid_Specific_Outpatient_Encounter_Existence), "No", Covid_Specific_Outpatient_Encounter_Existence),
           Covid_Specific_Outpatient_Consecutive_Encounters_Existence = ifelse(is.na(Covid_Specific_Outpatient_Consecutive_Encounters_Existence), "No", Covid_Specific_Outpatient_Consecutive_Encounters_Existence))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an outpatient encounter and status information"))
  
  rm(Outpatient, Covid_Outpatient, Covid_Outpatient_Status)
  
  # Compare Outpatient-Emergency metrics
  Emergency <- Encounters %>% filter(Inpatient_Outpatient == "Outpatient-Emergency")
  Output_columns <- Emergency %>% group_by(EMPI) %>% summarise(Covid_Related_with_Emergency_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Related_with_Emergency_Encounter_Info = ifelse(is.na(Covid_Related_with_Emergency_Encounter_Info), "No", Covid_Related_with_Emergency_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had an encounter that had a COVID diagnosis and at least one outpatient-emergency encounter"))
  
  Covid_Emergency <- Emergency %>% filter(contains_covid)
  Output_columns <- Covid_Emergency %>% group_by(EMPI) %>% summarise(Covid_Specific_Emergency_Encounter_Info = "Yes", .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Emergency_Encounter_Info = ifelse(is.na(Covid_Specific_Emergency_Encounter_Info), "No", Covid_Specific_Emergency_Encounter_Info))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an outpatient-emergency encounter"))
  
  Covid_Emergency_Status <- Covid_Emergency %>% filter(Encounter_Status == "Regular" |
                                                         Encounter_Status == "Inhouse" | 
                                                         Encounter_Status == "Contains Invalid DRGs") %>%
    arrange(EMPI, Admit_Date, desc(Discharge_Date))
  Covid_Emergency_Status <- distinct(Covid_Emergency_Status, EMPI, Admit_Date, .keep_all = TRUE)
  Covid_Emergency_Status <- distinct(Covid_Emergency_Status, EMPI, Discharge_Date, .keep_all = TRUE)
  Output_columns <- Covid_Emergency_Status %>% group_by(EMPI) %>%
    summarise(Covid_Specific_Emergency_Encounter_Existence = "Yes",
              Covid_Specific_Emergency_Encounter_Occurances = n(),
              Covid_Specific_Emergency_Consecutive_Encounters_Existence = ifelse(any(ifelse(any(LOS_Days > 0), TRUE, FALSE),
                                                                                     ifelse(Covid_Specific_Emergency_Encounter_Occurances > 1,
                                                                                            any(Admit_Date[2:length(Admit_Date)] - Discharge_Date[1:length(Discharge_Date) - 1] == 1),
                                                                                            FALSE)),
                                                                                 "Yes", "No"),
              Covid_Specific_Emergency_LOS_Total = sum(LOS_Days) + ifelse(Covid_Specific_Emergency_Encounter_Occurances > 1,
                                                                          sum(Admit_Date[2:length(Admit_Date)] - Discharge_Date[1:length(Discharge_Date) - 1] == 1),
                                                                          0),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Covid_Specific_Emergency_Encounter_Existence = ifelse(is.na(Covid_Specific_Emergency_Encounter_Existence), "No", Covid_Specific_Emergency_Encounter_Existence),
           Covid_Specific_Emergency_Consecutive_Encounters_Existence = ifelse(is.na(Covid_Specific_Emergency_Consecutive_Encounters_Existence), "No", Covid_Specific_Emergency_Consecutive_Encounters_Existence))
  loginfo(str_c(nrow(Output_columns), " subjects have had a COVID diagnosis during an outpatient-emergency encounter and status information"))
  
  rm(Emergency, Covid_Emergency, Covid_Emergency_Status)
  rm(Encounters, Output_columns)
  return(DF_to_fill)
}