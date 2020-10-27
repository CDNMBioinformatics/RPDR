require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_physical <- function(DF_to_fill = All_merged,
                             input_file_header = config$rpdr_file_header,
                             input_file_ending = config$rpdr_file_ending,
                             path_phy_abn = str_c(config$data_dir, "Phy_abnormalities/"),
                             output_file_ending = config$general_file_ending,
                             Return_BMI = TRUE,
                             Underweight_Normal = config$BMI_params$Underweight_Normal,
                             Normal_Overweight = config$BMI_params$Normal_Overweight,
                             Overweight_Obese = config$BMI_params$Overweight_Obese,
                             Return_Influenza = TRUE,
                             Return_Smoker = TRUE){
  loginfo("Processing Health History & Physical Findings data...")
  Phy <- data.table(fread(str_c(input_file_header, "Phy", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn)}
  logdebug(str_c("Note: All Health History & Physical Findings abnormalites can be found at ",
                 path_phy_abn))
  path_phy_abn = str_c(path_phy_abn, "Phy_")
  Phy <- Phy %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  
  if(Return_BMI){
    BMI <- Phy %>% filter(grepl("BMI", Concept_Name))
    # Underweight =                  [0 <= x < Underweight_Normal]
    # Normal      = [Underweight_Normal <= x < Normal_Overweight]
    # Overweight  =  [Normal_Overweight <= x < Overweight_Obese]
    # Obese       =   [Overweight_Obese <= x]
    # These values are based off of https://www.cancer.org/cancer/cancer-causes/diet-physical-activity/body-weight-and-cancer-risk/adult-bmi.html
    if(is.null(Underweight_Normal)){Underweight_Normal = 18.5}
    if(is.null(Normal_Overweight)){Normal_Overweight = 25}
    if(is.null(Overweight_Obese)){Overweight_Obese = 30}
    Phy_abn <- BMI[(duplicated(BMI)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(BMI), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_BMI_Duplicate_rows", output_file_ending))
      BMI <- BMI %>% unique()
    }
    rm(Phy_abn)
    BMI <- BMI %>% mutate(Result = as.numeric(Result))
    loginfo(str_c("Using ", Underweight_Normal, ", ", Normal_Overweight,", and ",
                  Overweight_Obese, " as thresholds between Underweight, Normal, Overweight, and Obese"))
    
    Output_Columns <- BMI %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI = "Yes",
                BMI_number_recorded = n(),
                BMI_average = mean(Result),
                BMI_min = min(Result, na.rm = TRUE),
                BMI_max = max(Result, na.rm = TRUE),
                BMI_median = median(Result, na.rm = TRUE),
                BMI_all_values = paste(Result, collapse = ";"),
                BMI_all_dates = paste(Date, collapse = ";"),
                BMI_probable_category = ifelse(BMI_median < Underweight_Normal, "Underweight",
                                               ifelse(BMI_median < Normal_Overweight, "Normal",
                                                      ifelse(BMI_median < Overweight_Obese,
                                                             "Overweight", "Obese"))),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>% mutate(BMI = ifelse(is.na(BMI), "No", BMI),
                                        BMI_number_recorded = ifelse(is.na(BMI_number_recorded),
                                                                     0, BMI_number_recorded))
    Underweight <- BMI %>% filter(Result < Underweight_Normal)
    Output_Columns <- Underweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_underweight = "Yes",
                BMI_underweight_number_recorded = n(),
                BMI_underweight_all_values = paste(Result, collapse = ";"),
                BMI_underweight_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_underweight = ifelse(is.na(BMI_underweight), "No", BMI_underweight),
             BMI_underweight_number_recorded = ifelse(is.na(BMI_underweight_number_recorded),
                                                      0, BMI_underweight_number_recorded))
    Normal <- BMI %>% filter(Result >= Underweight_Normal & Result < Normal_Overweight)
    Output_Columns <- Normal %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_normal = "Yes",
                BMI_normal_number_recorded = n(),
                BMI_normal_all_values = paste(Result, collapse = ";"),
                BMI_normal_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_normal = ifelse(is.na(BMI_normal), "No", BMI_normal),
             BMI_normal_number_recorded = ifelse(is.na(BMI_normal_number_recorded),
                                                 0, BMI_normal_number_recorded))
    Overweight <- BMI %>% filter(Result >= Normal_Overweight & Result < Overweight_Obese)
    Output_Columns <- Overweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_overweight = "Yes",
                BMI_overweight_number_recorded = n(),
                BMI_overweight_all_values = paste(Result, collapse = ";"),
                BMI_overweight_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_overweight = ifelse(is.na(BMI_overweight), "No", BMI_overweight),
             BMI_overweight_number_recorded = ifelse(is.na(BMI_overweight_number_recorded),
                                                     0, BMI_overweight_number_recorded))
    Obese <- BMI %>% filter(Result >= Overweight_Obese)
    Output_Columns <- Obese %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_obese = "Yes",
                BMI_obese_number_recorded = n(),
                BMI_obese_all_values = paste(Result, collapse = ";"),
                BMI_obese_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_obese = ifelse(is.na(BMI_obese), "No", BMI_obese),
             BMI_obese_number_recorded = ifelse(is.na(BMI_obese_number_recorded),
                                                0, BMI_obese_number_recorded))
    rm(Underweight, Normal, Overweight, Obese, BMI, Output_Columns)
  }
  if (Return_Influenza){
    Flu <- Phy %>% filter(grepl("Influenza", Concept_Name))
    Phy_abn <- Flu %>% filter(!grepl("^$|Done", Result))
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " out of ", nrow(Flu), " entries removed as records are of",
                    "declines, deferreds, and unavailables"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Flu_Invalid_information", output_file_ending))
      Flu <- Flu %>% filter(grepl("^$|Done", Result))
    }
    Phy_abn <- Flu[(duplicated(Flu)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(Flu), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Flu_Duplicate_rows", output_file_ending))
      Flu <- Flu %>% unique()
    }
    rm(Phy_abn)
    Output_Columns <- Flu %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(Flu_vaccined = "Yes",
                Flu_vaccine_count  = n(),
                Flu_vaccine_most_recent = last(Date),
                Flu_vaccine_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(Flu_vaccined = ifelse(is.na(Flu_vaccined), "No", Flu_vaccined),
             Flu_vaccine_count = ifelse(is.na(Flu_vaccine_count), 0, Flu_vaccine_count))
    rm(Output_Columns, Flu)
  }
  if(Return_Smoker){
    Smoker = Phy %>% filter(Concept_Name %in% c("Smoking Quit Date",
                                                "Smoking Start Date",
                                                "Smoking Tobacco Use-Current Every Day Smoker",
                                                "Smoking Tobacco Use-Current Some Day Smoker",
                                                "Smoking Tobacco Use-Former Smoker",
                                                "Smoking Tobacco Use-Heavy Tobacco Smoker",
                                                "Smoking Tobacco Use-Light Tobacco Smoker",
                                                "Smoking Tobacco Use-Smoker, Current Status Unknown"))
    Phy_abn <- Smoker[(duplicated(Smoker)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(Smoker), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Smoker_Duplicate_rows", output_file_ending))
      Smoker <- Smoker %>% unique()
    }
    rm(Phy_abn)
    Output_Columns <- Smoker %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(Smoker_Former_or_Current = "Yes",
                Smoker_Former_or_Current_First_Date = first(Date),
                Smoker_Former_or_Current_Most_Recent_Date = last(Date),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(Smoker_Former_or_Current = ifelse(is.na(Smoker_Former_or_Current), "No", Smoker_Former_or_Current))
    rm(Output_Columns, Smoker)
  }
  return(DF_to_fill)
}
process_physical_date <- function(DF_to_fill = All_merged,
                                  input_file_header = config$rpdr_file_header,
                                  input_file_ending = config$rpdr_file_ending,
                                  path_phy_abn = str_c(config$data_dir, "Phy_abnormalities/"),
                                  output_file_ending = config$general_file_ending,
                                  Return_BMI = TRUE,
                                  Underweight_Normal = config$BMI_params$Underweight_Normal,
                                  Normal_Overweight = config$BMI_params$Normal_Overweight,
                                  Overweight_Obese = config$BMI_params$Overweight_Obese,
                                  Return_Influenza = TRUE,
                                  Return_Smoker = TRUE,
                                  Date_Column,
                                  restrict_to_before = FALSE){
  loginfo("Processing Health History & Physical Findings data...")
  Phy <- fread(str_c(input_file_header, "Phy", input_file_ending)) %>% arrange(EMPI)
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn)}
  logdebug(str_c("Note: All Health History & Physical Findings abnormalites can be found at ",
                 path_phy_abn))
  path_phy_abn <- str_c(path_phy_abn, ifelse(restrict_to_before,
                                             "Phy_PreDate_",
                                             "Phy_PostDate_"))
  if (!missing(Date_Column)){
    Relevant_rows <- NULL
    Original_Columns <- names(DF_to_fill)
  }
  loginfo("Restricting prescription by cutoff...")
  # Restrict ids to only the ones in the reduced id list
  Phy <- Phy %>% filter(EMPI %in% DF_to_fill$EMPI)
  Phy <- Phy %>% arrange(EMPI, Date)
  if (!sum(str_count(DF_to_fill %>% pull(Date_Column)) == 10) == nrow(DF_to_fill)){
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, Date_Column) %>%
      rename(Cutoff_Date_Time = Date_Column) %>%
      extract(Cutoff_Date_Time, c("Cutoff_Date", "Cutoff_Time"),
              regex = "(\\d{4}-\\d{2}-\\d{2}) (\\d{2}:\\d{2})", remove = FALSE) %>%
      mutate(Cutoff_Date = ifelse(is.na(Cutoff_Time), Cutoff_Date_Time, Cutoff_Date)) %>% select(EMPI, Cutoff_Date)
  } else {
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, Date_Column) %>%
      rename(Cutoff_Date = Date_Column)
  }
  # Restrict Phy to only before or only after
  Phy <- left_join(Phy, EMPI_Date_Limit, by = "EMPI")
  time_str = ifelse(restrict_to_before, "before", "after")
  loginfo(str_c("Restricting Phy to ", time_str, " cutoff dates... "))
  
  Original_row_length <- nrow(Phy)
  if(restrict_to_before){
    Phy <- Phy %>% filter(Date <= Cutoff_Date)
  } else {
    Phy <- Phy %>% filter(Date >= Cutoff_Date)
  }
  New_row_length <- nrow(Phy)
  logdebug(str_c(Original_row_length, " Phy reduced to ", New_row_length, " Phy"))
  rm(time_str, Original_row_length, New_row_length)
  
  if(Return_BMI){
    BMI <- Phy %>% filter(grepl("BMI", Concept_Name))
    # Underweight =                  [0 <= x < Underweight_Normal]
    # Normal      = [Underweight_Normal <= x < Normal_Overweight]
    # Overweight  =  [Normal_Overweight <= x < Overweight_Obese]
    # Obese       =   [Overweight_Obese <= x]
    # These values are based off of https://www.cancer.org/cancer/cancer-causes/diet-physical-activity/body-weight-and-cancer-risk/adult-bmi.html
    if(is.null(Underweight_Normal)){Underweight_Normal = 18.5}
    if(is.null(Normal_Overweight)){Normal_Overweight = 25}
    if(is.null(Overweight_Obese)){Overweight_Obese = 30}
    Phy_abn <- BMI[(duplicated(BMI)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(BMI), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_BMI_Duplicate_rows", output_file_ending))
      BMI <- BMI %>% unique()
    }
    rm(Phy_abn)
    BMI <- BMI %>% mutate(Result = as.numeric(Result))
    loginfo(str_c("Using ", Underweight_Normal, ", ", Normal_Overweight,", and ",
                  Overweight_Obese, " as thresholds between Underweight, Normal, Overweight, and Obese"))
    
    Output_Columns <- BMI %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI = "Yes",
                BMI_number_recorded = n(),
                BMI_average = mean(Result),
                BMI_min = min(Result, na.rm = TRUE),
                BMI_max = max(Result, na.rm = TRUE),
                BMI_median = median(Result, na.rm = TRUE),
                BMI_all_values = paste(Result, collapse = ";"),
                BMI_all_dates = paste(Date, collapse = ";"),
                BMI_probable_category = ifelse(BMI_median < Underweight_Normal, "Underweight",
                                               ifelse(BMI_median < Normal_Overweight, "Normal",
                                                      ifelse(BMI_median < Overweight_Obese,
                                                             "Overweight", "Obese"))),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>% mutate(BMI = ifelse(is.na(BMI), "No", BMI),
                                        BMI_number_recorded = ifelse(is.na(BMI_number_recorded),
                                                                     0, BMI_number_recorded))
    Underweight <- BMI %>% filter(Result < Underweight_Normal)
    Output_Columns <- Underweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_underweight = "Yes",
                BMI_underweight_number_recorded = n(),
                BMI_underweight_all_values = paste(Result, collapse = ";"),
                BMI_underweight_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_underweight = ifelse(is.na(BMI_underweight), "No", BMI_underweight),
             BMI_underweight_number_recorded = ifelse(is.na(BMI_underweight_number_recorded),
                                                      0, BMI_underweight_number_recorded))
    Normal <- BMI %>% filter(Result >= Underweight_Normal & Result < Normal_Overweight)
    Output_Columns <- Normal %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_normal = "Yes",
                BMI_normal_number_recorded = n(),
                BMI_normal_all_values = paste(Result, collapse = ";"),
                BMI_normal_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_normal = ifelse(is.na(BMI_normal), "No", BMI_normal),
             BMI_normal_number_recorded = ifelse(is.na(BMI_normal_number_recorded),
                                                 0, BMI_normal_number_recorded))
    Overweight <- BMI %>% filter(Result >= Normal_Overweight & Result < Overweight_Obese)
    Output_Columns <- Overweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_overweight = "Yes",
                BMI_overweight_number_recorded = n(),
                BMI_overweight_all_values = paste(Result, collapse = ";"),
                BMI_overweight_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_overweight = ifelse(is.na(BMI_overweight), "No", BMI_overweight),
             BMI_overweight_number_recorded = ifelse(is.na(BMI_overweight_number_recorded),
                                                     0, BMI_overweight_number_recorded))
    Obese <- BMI %>% filter(Result >= Overweight_Obese)
    Output_Columns <- Obese %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_obese = "Yes",
                BMI_obese_number_recorded = n(),
                BMI_obese_all_values = paste(Result, collapse = ";"),
                BMI_obese_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_obese = ifelse(is.na(BMI_obese), "No", BMI_obese),
             BMI_obese_number_recorded = ifelse(is.na(BMI_obese_number_recorded),
                                                0, BMI_obese_number_recorded))
    rm(Underweight, Normal, Overweight, Obese, BMI, Output_Columns)
  }
  if (Return_Influenza){
    Flu <- Phy %>% filter(grepl("Influenza", Concept_Name))
    Phy_abn <- Flu %>% filter(!grepl("^$|Done", Result))
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " out of ", nrow(Flu), " entries removed as records are of",
                    "declines, deferreds, and unavailables"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Flu_Invalid_information", output_file_ending))
      Flu <- Flu %>% filter(grepl("^$|Done", Result))
    }
    Phy_abn <- Flu[(duplicated(Flu)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(Flu), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Flu_Duplicate_rows", output_file_ending))
      Flu <- Flu %>% unique()
    }
    rm(Phy_abn)
    Output_Columns <- Flu %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(Flu_vaccined = "Yes",
                Flu_vaccine_count  = n(),
                Flu_vaccine_most_recent = last(Date),
                Flu_vaccine_all_dates = paste(Date, collapse = ";"),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(Flu_vaccined = ifelse(is.na(Flu_vaccined), "No", Flu_vaccined),
             Flu_vaccine_count = ifelse(is.na(Flu_vaccine_count), 0, Flu_vaccine_count))
    rm(Output_Columns, Flu)
  }
  if(Return_Smoker){
    Smoker = Phy %>% filter(Concept_Name %in% c("Smoking Quit Date",
                                                "Smoking Start Date",
                                                "Smoking Tobacco Use-Current Every Day Smoker",
                                                "Smoking Tobacco Use-Current Some Day Smoker",
                                                "Smoking Tobacco Use-Former Smoker",
                                                "Smoking Tobacco Use-Heavy Tobacco Smoker",
                                                "Smoking Tobacco Use-Light Tobacco Smoker",
                                                "Smoking Tobacco Use-Smoker, Current Status Unknown"))
    Phy_abn <- Smoker[(duplicated(Smoker)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(Smoker), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_Smoker_Duplicate_rows", output_file_ending))
      Smoker <- Smoker %>% unique()
    }
    rm(Phy_abn)
    Output_Columns <- Smoker %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(Smoker_Former_or_Current = "Yes",
                Smoker_Former_or_Current_First_Date = first(Date),
                Smoker_Former_or_Current_Most_Recent_Date = last(Date),
                .groups = "drop")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(Smoker_Former_or_Current = ifelse(is.na(Smoker_Former_or_Current), "No", Smoker_Former_or_Current))
    rm(Output_Columns, Smoker)
  }
  return(DF_to_fill)
}