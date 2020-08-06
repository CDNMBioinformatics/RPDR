require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process7_physical <- function(DF_to_fill = All_merged,
                              input_file_header = config$rpdr_file_header,
                              input_file_ending = config$rpdr_file_ending,
                              path_phy_abn = str_c(config$data_dir, "Phy_abnormalities/"),
                              output_file_ending = config$general_file_ending,
                              Return_BMI = TRUE,
                              Underweight_Normal = 18.5,
                              Normal_Overweight = 24.9,
                              Overweight_Obese = 30,
                              Return_Influenza = TRUE){
  loginfo("Processing Health History & Physical Findings data...")
  Phy <- data.table(fread(str_c(input_file_header, "Phy", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn)}
  logdebug(str_c("Note: All Health History & Physical Findings abnormalites can be found at ",
                 path_phy_abn))
  if(Return_BMI){
    BMI <- Phy %>% filter(grepl("BMI", Concept_Name))
    
    Phy_abn <- BMI[(duplicated(BMI)),]
    if (nrow(Phy_abn) > 0){
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(BMI), " removed"))
      fwrite(Phy_abn, str_c(path_phy_abn, "Abnormality_BMI_Duplicate_rows", output_file_ending))
      BMI <- BMI %>% unique()
    }
    rm(Phy_abn)
    BMI <- BMI %>% mutate(Result = as.numeric(Result),
                          Date = mdy(Date))
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
                                                             "Overweight", "Obese"))))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
    DF_to_fill <- DF_to_fill %>% mutate(BMI = ifelse(is.na(BMI), "No", BMI),
                                        BMI_number_recorded = ifelse(is.na(BMI_number_recorded),
                                                                     0, BMI_number_recorded))
    Underweight <- BMI %>% filter(Result < Underweight_Normal)
    Output_Columns <- Underweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_underweight = "Yes",
                BMI_underweight_number_recorded = n(),
                BMI_underweight_all_values = paste(Result, collapse = ";"),
                BMI_underweight_all_dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_underweight = ifelse(is.na(BMI_underweight), "No", BMI_underweight),
             BMI_underweight_number_recorded = ifelse(is.na(BMI_underweight_number_recorded),
                                                      0, BMI_underweight_number_recorded))
    Normal <- BMI %>% filter(Result >= Underweight_Normal & Result < Normal_Overweight)
    Output_Columns <- Normal %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_normal = "Yes",
                BMI_normal_number_recorded = n(),
                BMI_normal_all_values = paste(Result, collapse = ";"),
                BMI_normal_all_dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_normal = ifelse(is.na(BMI_normal), "No", BMI_normal),
             BMI_normal_number_recorded = ifelse(is.na(BMI_normal_number_recorded),
                                                 0, BMI_normal_number_recorded))
    Overweight <- BMI %>% filter(Result >= Normal_Overweight & Result < Overweight_Obese)
    Output_Columns <- Overweight %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_overweight = "Yes",
                BMI_overweight_number_recorded = n(),
                BMI_overweight_all_values = paste(Result, collapse = ";"),
                BMI_overweight_all_dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(BMI_overweight = ifelse(is.na(BMI_overweight), "No", BMI_overweight),
             BMI_overweight_number_recorded = ifelse(is.na(BMI_overweight_number_recorded),
                                                     0, BMI_overweight_number_recorded))
    Obese <- BMI %>% filter(Result >= Overweight_Obese)
    Output_Columns <- Obese %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(BMI_obese = "Yes",
                BMI_obese_number_recorded = n(),
                BMI_obese_all_values = paste(Result, collapse = ";"),
                BMI_obese_all_dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
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
    Flu <- Flu %>% mutate(Date = mdy(Date))
    Output_Columns <- Flu %>% group_by(EMPI) %>% arrange(Date) %>%
      summarise(Flu_vaccined = "Yes",
                Flu_vaccine_count  = n(),
                Flu_vaccine_most_recent = last(Date),
                Flu_vaccine_all_dates = paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns)
    DF_to_fill <- DF_to_fill %>%
      mutate(Flu_vaccined = ifelse(is.na(Flu_vaccined), "No", Flu_vaccined),
             Flu_vaccine_count = ifelse(is.na(Flu_vaccine_count), 0, Flu_vaccine_count))
    rm(Output_Columns, Flu)
  }
  return(DF_to_fill)
}