require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_procedures <- function(DF_to_fill = All_merged,
                               input_file_header = config$rpdr_file_header,
                               input_file_ending = config$rpdr_file_ending,
                               path_prc_abn = str_c(config$data_dir, "Prc_abnormalities/"),
                               output_file_ending = config$general_file_ending){
  loginfo("Processing Procedures data...")
  Procedures <- data.table(fread(str_c(input_file_header, "Prc", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_prc_abn)) {dir.create(path_prc_abn)}
  logdebug(str_c("Note: All procedure abnormalites can be found at ", path_prc_abn))
  Procedures <- Procedures %>% mutate(Date = as.Date(Date, "%m/%d/%Y")) %>% arrange(EMPI, Date)
  
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
  Procedures <- left_join(Procedures, Dates_of_Interest)
  
  vent <- Procedures %>% filter(grepl("[Vv]entilat", Procedure_Name))
  
  Mechanical <- Procedures %>% filter(grepl("mechanical ventilation", Procedure_Name))
  Output_columns <- Mechanical %>% group_by(EMPI) %>%
    summarise(Any_Mechanical_Ventilation = "Yes",
              Any_Mechanical_Ventilation_Count = n(),
              Any_Mechanical_Ventilation_Dates = paste(Date, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Any_Mechanical_Ventilation = ifelse(is.na(Any_Mechanical_Ventilation),
                                               "No", Any_Mechanical_Ventilation))
  
  Output_columns <- Mechanical %>% filter(grepl("more", Procedure_Name)) %>% group_by(EMPI) %>%
    summarise(Any_Mechanical_Ventilation_GTEq_96_Hours = "Yes",
              Any_Mechanical_Ventilation_GTEq_96_Hours_Count = n(),
              Any_Mechanical_Ventilation_GTEq_96_Hours_Dates = paste(Date, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Any_Mechanical_Ventilation = ifelse(is.na(Any_Mechanical_Ventilation_GTEq_96_Hours),
                                               "No", Any_Mechanical_Ventilation_GTEq_96_Hours))
  
  Output_columns <- Mechanical %>% filter(grepl("less", Procedure_Name)) %>% group_by(EMPI) %>%
    summarise(Any_Mechanical_Ventilation_LT_96_Hours = "Yes",
              Any_Mechanical_Ventilation_LT_96_Hours_Count = n(),
              Any_Mechanical_Ventilation_LT_96_Hours_Dates = paste(Date, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Any_Mechanical_Ventilation = ifelse(is.na(Any_Mechanical_Ventilation_LT_96_Hours),
                                               "No", Any_Mechanical_Ventilation_LT_96_Hours))
  
  Output_columns <- Mechanical %>% filter(grepl("unspecified", Procedure_Name)) %>% group_by(EMPI) %>%
    summarise(Any_Mechanical_Ventilation_Unknown_Hours = "Yes",
              Any_Mechanical_Ventilation_Unknown_Hours_Count = n(),
              Any_Mechanical_Ventilation_Unknown_Hours_Dates = paste(Date, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Any_Mechanical_Ventilation = ifelse(is.na(Any_Mechanical_Ventilation_Unknown_Hours),
                                               "No", Any_Mechanical_Ventilation_Unknown_Hours))
  
  Mechanical <- Mechanical %>% filter(Date >= Lower_Bound_Date)
  if (nrow(Mechanical)){
    Output_columns <- Mechanical %>% group_by(EMPI) %>%
      summarise(Mechanical_Ventilation_After_Covid = "Yes",
                Mechanical_Ventilation_Count = n(),
                Mechanical_Ventilation_After_Covid_Dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(Mechanical_Ventilation_After_Covid = ifelse(is.na(Mechanical_Ventilation_After_Covid),
                                                         "No", Mechanical_Ventilation_After_Covid))
    
    Output_columns <- Mechanical %>% filter(grepl("more", Procedure_Name)) %>% group_by(EMPI) %>%
      summarise(Mechanical_Ventilation_After_Covid_GTEq_96_Hours = "Yes",
                Mechanical_Ventilation_After_Covid_GTEq_96_Hours_Count = n(),
                Mechanical_Ventilation_After_Covid_GTEq_96_Hours_Dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(Mechanical_Ventilation_After_Covid = ifelse(is.na(Mechanical_Ventilation_After_Covid_GTEq_96_Hours),
                                                         "No", Mechanical_Ventilation_After_Covid_GTEq_96_Hours))
    
    Output_columns <- Mechanical %>% filter(grepl("less", Procedure_Name)) %>% group_by(EMPI) %>%
      summarise(Mechanical_Ventilation_After_Covid_LT_96_Hours = "Yes",
                Mechanical_Ventilation_After_Covid_LT_96_Hours_Count = n(),
                Mechanical_Ventilation_After_Covid_LT_96_Hours_Dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(Mechanical_Ventilation_After_Covid = ifelse(is.na(Mechanical_Ventilation_After_Covid_LT_96_Hours),
                                                         "No", Mechanical_Ventilation_After_Covid_LT_96_Hours))
    
    Output_columns <- Mechanical %>% filter(grepl("unspecified", Procedure_Name)) %>% group_by(EMPI) %>%
      summarise(Mechanical_Ventilation_After_Covid_Unknown_Hours = "Yes",
                Mechanical_Ventilation_After_Covid_Unknown_Hours_Count = n(),
                Mechanical_Ventilation_After_Covid_Unknown_Hours_Dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(Mechanical_Ventilation_After_Covid = ifelse(is.na(Mechanical_Ventilation_After_Covid_Unknown_Hours),
                                                         "No", Mechanical_Ventilation_After_Covid_Unknown_Hours))
  } else {
    loginfo("No Mechanical Ventilation usage found after covid")
  }
  return(DF_to_fill)
}