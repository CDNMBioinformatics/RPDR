require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

####################
# Helper Functions #
####################
HealthHistoryHelper_BMIOutputs <- function(BMI_DF,
                                           Result_DF,
                                           BMI_cat_name){
  BMI_Cat <- str_c("BMI_", BMI_cat_name)
  BMI_Cat_Num <- str_c(BMI_Cat, "_number_recorded")
  BMI_Cat_Val <- str_c(BMI_Cat, "_all_values")
  BMI_Cat_Dat <- str_c(BMI_Cat, "_all_dates")
  Output_Columns <- BMI_DF %>%	
    group_by(EMPI) %>% arrange(Date) %>%	
    summarise(!!as.symbol(BMI_Cat) := "Yes",	
              !!as.symbol(BMI_Cat_Num) := n(),	
              !!as.symbol(BMI_Cat_Val) := paste(Result, collapse = ";"),	
              !!as.symbol(BMI_Cat_Dat) := paste(Date, collapse = ";"),	
              .groups = "drop")	
  Result_DF <- left_join(Result_DF, Output_Columns, by = "EMPI")	
  Result_DF <- Result_DF %>%	
    mutate(!!as.symbol(BMI_Cat) := ifelse(is.na(!!as.symbol(BMI_Cat)), "No", !!as.symbol(BMI_Cat)),	
           !!as.symbol(BMI_Cat_Num) := ifelse(is.na(!!as.symbol(BMI_Cat_Num)), 0,	!!as.symbol(BMI_Cat_Num)))	
  return(Result_DF)
}

HealthHisotryHelper_BPOutputs <- function(BP_DF,
                                          Result_DF,
                                          BP_cat_name){
  BP_name <- str_c("BP_", gsub(":", "", gsub(" ", "_", BP_cat_name)))
  BP_number <- str_c(BP_name, "_number_recorded")
  Output_Columns <- BP_DF %>%	
    group_by(EMPI) %>%	
    filter(Category == BP_cat_name) %>%	
    summarise(!!as.symbol(BP_name) := "Yes",
              !!as.symbol(BP_number) := n())	
  Result_DF <- left_join(Result_DF, Output_Columns, by = "EMPI")	
  Result_DF <- Result_DF %>% 
    mutate(!!as.symbol(BP_name) := ifelse(is.na(!!as.symbol(BP_name)), "No", !!as.symbol(BP_name)),	
           !!as.symbol(BP_number) := ifelse(is.na(!!as.symbol(BP_number)), 0, !!as.symbol(BP_number)))
  return(Result_DF)
}

HealthHistoryHelper_CoreFunctionality <- function(Physical_DF,
                                                  Main_DF,
                                                  Path_Abnormality,
                                                  File_Suffix,
                                                  Gen_BMI,
                                                  Und_Norm,
                                                  Norm_Over,
                                                  Over_Ob,
                                                  Gen_Flu,
                                                  Gen_Smoke,
                                                  Gen_BP){
  if (Gen_BMI){	
    BMI <- Physical_DF %>% filter(grepl("BMI", Concept_Name))	
    # Underweight =         [0 <= x < Und_Norm]	
    # Normal      =  [Und_Norm <= x < Norm_Over]	
    # Overweight  = [Norm_Over <= x < Over_Ob]	
    # Obese       =   [Over_Ob <= x]	
    # These values are based off of https://www.cancer.org/cancer/cancer-causes/diet-physical-activity/body-weight-and-cancer-risk/adult-bmi.html	
    if(is.null(Und_Norm)){Und_Norm = 18.5}	
    if(is.null(Norm_Over)){Norm_Over = 25}	
    if(is.null(Over_Ob)){Over_Ob = 30}
    Phy_abn <- BMI[(duplicated(BMI)),]	
    if (nrow(Phy_abn) > 0){	
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(BMI), " removed"))	
      fwrite(Phy_abn, str_c(Path_Abnormality, "Abnormality_BMI_Duplicate_rows", File_Suffix))	
      BMI <- BMI %>% unique()	
    }	
    rm(Phy_abn)	
    BMI <- BMI %>% mutate(Result = as.numeric(Result))	
    loginfo(str_c("Using ", Und_Norm, ", ", Norm_Over,", and ",	
                  Over_Ob, " as thresholds between Underweight, Normal, Overweight, and Obese"))	
    
    Output_Columns <- BMI %>% group_by(EMPI) %>% arrange(Date) %>%	
      summarise(BMI = "Yes",	
                BMI_number_recorded = n(),	
                BMI_average = mean(Result),	
                BMI_min = min(Result, na.rm = TRUE),	
                BMI_max = max(Result, na.rm = TRUE),	
                BMI_median = median(Result, na.rm = TRUE),	
                BMI_all_values = paste(Result, collapse = ";"),	
                BMI_all_dates = paste(Date, collapse = ";"),	
                BMI_probable_category = case_when(BMI_median < Und_Norm ~ "Underweight",	
                                                  BMI_median < Norm_Over ~ "Normal",	
                                                  BMI_median < Over_Ob ~ "Overweight",	
                                                  TRUE ~ "Obese"),	
                .groups = "drop")	
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")	
    Main_DF <- Main_DF %>%	
      mutate(BMI = ifelse(is.na(BMI), "No", BMI),	
             BMI_number_recorded = ifelse(is.na(BMI_number_recorded),	0, BMI_number_recorded))	
    Underweight <- BMI %>% filter(Result < Und_Norm)
    Main_DF <- HealthHistoryHelper_BMIOutputs(BMI_DF = Underweight,
                                              Result_DF = Main_DF,
                                              BMI_cat_name = "underweight")
    Normal <- BMI %>% filter(Result >= Und_Norm & Result < Norm_Over)
    Main_DF <- HealthHistoryHelper_BMIOutputs(BMI_DF = Normal,
                                              Result_DF = Main_DF,
                                              BMI_cat_name = "normal")
    Overweight <- BMI %>% filter(Result >= Norm_Over & Result < Over_Ob)
    Main_DF <- HealthHistoryHelper_BMIOutputs(BMI_DF = Overweight,
                                              Result_DF = Main_DF,
                                              BMI_cat_name = "overweight")
    Obese <- BMI %>% filter(Result >= Over_Ob)
    Main_DF <- HealthHistoryHelper_BMIOutputs(BMI_DF = Obese,
                                              Result_DF = Main_DF,
                                              BMI_cat_name = "obese")
    rm(BMI, Underweight, Normal, Overweight, Obese, Output_Columns)	
  }
  if (Gen_Flu){	
    Flu <- Physical_DF %>% filter(grepl("Influenza", Concept_Name))	
    Phy_abn <- Flu %>% filter(!grepl("^$|Done", Result))	
    if (nrow(Phy_abn) > 0){	
      logwarn(str_c(nrow(Phy_abn), " out of ", nrow(Flu), " entries removed as records are of",	
                    "declines, deferreds, and unavailables"))	
      fwrite(Phy_abn, str_c(Path_Abnormality, "Abnormality_Flu_Invalid_information", File_Suffix))	
      Flu <- Flu %>% filter(grepl("^$|Done", Result))	
    }
    rm(Phy_abn)
    Phy_abn <- Flu[(duplicated(Flu)),]	
    if (nrow(Phy_abn) > 0){	
      logwarn(str_c(nrow(Phy_abn), " completely duplicated row(s) out of ", nrow(Flu), " removed"))	
      fwrite(Phy_abn, str_c(Path_Abnormality, "Abnormality_Flu_Duplicate_rows", File_Suffix))	
      Flu <- Flu %>% unique()	
    }	
    rm(Phy_abn)	
    Output_Columns <- Flu %>%	
      group_by(EMPI) %>% arrange(Date) %>%	
      summarise(Flu_vaccined = "Yes",	
                Flu_vaccine_count  = n(),	
                Flu_vaccine_most_recent = last(Date),	
                Flu_vaccine_all_dates = paste(Date, collapse = ";"),	
                .groups = "drop")	
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")	
    Main_DF <- Main_DF %>%	
      mutate(Flu_vaccined = ifelse(is.na(Flu_vaccined), "No", Flu_vaccined),	
             Flu_vaccine_count = ifelse(is.na(Flu_vaccine_count), 0, Flu_vaccine_count))	
    rm(Output_Columns, Flu)	
  }
  if (Gen_Smoke){	
    Smoker = Physical_DF %>% 
      filter(Concept_Name %in% c("Smoking Quit Date",	
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
      fwrite(Phy_abn, str_c(Path_Abnormality, "Abnormality_Smoker_Duplicate_rows", File_Suffix))	
      Smoker <- Smoker %>% unique()	
    }	
    rm(Phy_abn)	
    Output_Columns <- Smoker %>% group_by(EMPI) %>% arrange(Date) %>%	
      summarise(Smoker_Former_or_Current = "Yes",	
                Smoker_Former_or_Current_First_Date = first(Date),	
                Smoker_Former_or_Current_Most_Recent_Date = last(Date),	
                .groups = "drop")	
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")	
    Main_DF <- Main_DF %>%	
      mutate(Smoker_Former_or_Current = ifelse(is.na(Smoker_Former_or_Current),	
                                               "No",	
                                               Smoker_Former_or_Current))	
    rm(Output_Columns, Smoker)	
  }
  if (Gen_BP){
    #####################
    # Shorthand strings #
    #####################
    BP_List = list()
    BP_List$N = "Normal"
    BP_List$E = "Elevated"
    BP_List$HBP1 = "High Blood Pressure: Stage 1"
    BP_List$HBP2 = "High Blood Pressure: Stage 2"
    BP_List$HC = "Hypertensive Crisis"
    BP <- Physical_DF %>% filter(grepl("Blood [Pp]", Concept_Name))	
    # readings <- BP %>%	
    #   group_by(EMPI, Date) %>%	
    #   summarise(count = n(), results = paste(Result, collapse = ";"))	
    # readings_pure <- BP %>%	
    #   filter(grepl("^\\d", Result)) %>%	
    #   group_by(EMPI, Date) %>%	
    #   summarise(count = n(), results = paste(Result, collapse = ";"))	
    # BP %>% filter(grepl("^\\d+[^/]\\d+$", Result))	
    # Physical_DF %>% filter(grepl("Dias|Syst", Concept_Name)) %>% group_by(Concept_Name, Code) %>% summarise(n())	
    # There is a log of logic in cleaning this data that came with discussions from MDs and RNs about what errors likely came from (most are entry errors)
    BP <- Physical_DF %>%	
      filter(grepl("Blood [Pp]|Systolic/Diastolic", Concept_Name),	
             grepl("^\\d+/\\d+$", Result)) %>%	
      separate(Result, c("Systolic", "Diastolic"), sep = "/", remove = FALSE) %>%	
      mutate(Systolic = as.numeric(Systolic),	
             Diastolic = as.numeric(Diastolic)) %>%	
      mutate(Pulse_Pressure = Systolic - Diastolic,	
             PP_Check = round(Pulse_Pressure/Systolic, 2)) %>% # This should be around 1/3	
      # if BP is 5digits/2 digits it should actual first3digits/last2digits and the /2 digits is a different metric that got clumped in	
      mutate(Diastolic = ifelse(floor(log10(Systolic)) + 1 > 4,	
                                Systolic %% 100,	
                                Diastolic),	
             Systolic = ifelse(floor(log10(Systolic)) + 1 > 4,	
                               floor(Systolic/100),	
                               Systolic)) %>%	
      # A Diastolic less that 10 is not a thing (You are very likely dead at that point). Probably lost a zero along the way.	
      mutate(Diastolic = ifelse(floor(log10(Diastolic)) + 1 == 1,	
                                Diastolic*10,	
                                Diastolic)) %>%	
      # A Systolic more than 300 is not a thing (you are definitely dead at that point). Probably a weird truncation.	
      mutate(Systolic = ifelse(floor(log10(Systolic)) + 1 > 3,	
                               floor(Systolic/10),	
                               Systolic)) %>%	
      # Systolic should be greater than Diastolic	
      mutate(Diastolic = ifelse(Pulse_Pressure < 0 & Diastolic >= 200,	
                                Diastolic/10,	
                                Diastolic),	
             Systolic = ifelse(Pulse_Pressure < 0 & Systolic <= 50,	
                               Systolic*10,	
                               Systolic)) %>%	
      mutate(Pulse_Pressure = Systolic - Diastolic,	
             PP_Check = round(Pulse_Pressure/Systolic, 2)) %>% # This should be around 1/3	
      # If at this point there's still values that have a negative Pulse_Pressure, it's an entry error that's not easily fixable	
      filter(Pulse_Pressure > 0)
    # Categorizations based on https://www.heart.org/en/health-topics/high-blood-pressure/understanding-blood-pressure-readings
    BP <- BP %>%	
      mutate(Category = case_when(Systolic < 120 & Diastolic < 80 ~ BP_List$N,	
                                  Systolic < 130 & Diastolic < 80 ~ BP_List$E,	
                                  (Systolic >= 130 & Systolic < 140) |	
                                    (Diastolic >= 80 & Diastolic < 90) ~ BP_List$HBP1,	
                                  (Systolic >= 140 & Systolic < 180) |	
                                    (Diastolic >= 90 & Diastolic < 120) ~ BP_List$HBP2,	
                                  Systolic >=180 | Diastolic >= 120 ~ BP_List$HC,	
                                  TRUE ~ "Other"))	
    # Output_Columns <- BP %>%	
    #   group_by(EMPI, Category) %>% arrange(Date) %>%	
    #   summarise(n())	
    Output_Columns <- BP %>%	
      group_by(EMPI) %>% arrange(Date) %>%	
      summarise(BP = "Yes",	
                BP_number_recorded = n())	
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")	
    Main_DF <- Main_DF %>%	
      mutate(BP = ifelse(is.na(BP), "No", BP),	
             BP_number_recorded = ifelse(is.na(BP_number_recorded), 0, BP_number_recorded))
    for (BP_cat in unlist(BP_List)){
      Main_DF <- HealthHisotryHelper_BPOutputs(BP_DF = BP,
                                               Result_DF = Main_DF,
                                               BP_cat_name = BP_cat)
    }
    Output_Columns <- Main_DF %>%	
      filter(BP == "Yes") %>%	
      mutate(BP_largest_count = pmax(BP_Normal_number_recorded,	
                                     BP_Elevated_number_recorded,	
                                     BP_High_Blood_Pressure_Stage_1_number_recorded,	
                                     BP_High_Blood_Pressure_Stage_2_number_recorded,	
                                     BP_Hypertensive_Crisis_number_recorded,	
                                     na.rm = TRUE),	
             BP_probable_category = case_when(BP_largest_count == BP_Normal_number_recorded ~ BP_List$N,	
                                              BP_largest_count == BP_Elevated_number_recorded ~ BP_List$E,	
                                              BP_largest_count == BP_High_Blood_Pressure_Stage_1_number_recorded ~ BP_List$HBP1,	
                                              BP_largest_count == BP_High_Blood_Pressure_Stage_2_number_recorded ~ BP_List$HBP2,	
                                              TRUE ~ BP_List$HC))	
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")	
    rm(BP, Output_Columns)
  }
  return(Main_DF)
}

##################
# Main Functions #
##################
#' Read in file with diagnosis information and add relevant columns to existing data frame
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_phy_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Phy_abnormalities/)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param Return_BMI Do you want information relating to BMI? (default = TRUE)
#' @param Underweight_Normal The value that splits BMI categorization underweight and normal (default = config$BMI_params$Underweight_Normal)
#' @param Normal_Overweight The value that splits BMI categorization normal and overweight (default = config$BMI_params$Normal_Overweight)
#' @param Overweight_Obese The value that splits BMI categorization overweight and obese (default = config$BMI_params$Overweight_Obese)
#' @param Return_Influenza Do you want information relating to influenza? (default = TRUE)
#' @param Return_Smoker Do you want information relating to smoking history? (default = TRUE)
#' @param Return_Blood_Pressure Do you want information relating to blood pressure? (default = FALSE)
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_physical()
#' process_physical(Return_BMI = FALSE, Return_Blood_Pressure = TRUE)
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
                             Return_Smoker = TRUE,
                             Return_Blood_Pressure = FALSE){
  loginfo("Processing Health History & Physical Findings...")
  ###############################
  # Setup and Information Steps #
  ###############################
  Phy_file_name = str_c(input_file_header, "Phy", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Phy_file_name))
  
  loginfo("Generate columns relating to BMI: ", Return_BMI)
  if (Return_BMI){
    loginfo(str_c("Underweight: [BMI < ", Underweight_Normal, "]"))
    loginfo(str_c("Normal: [", Underweight_Normal, " <= BMI < ", Normal_Overweight, "]"))
    loginfo(str_c("Normal: [", Normal_Overweight, " <= BMI < ", Overweight_Obese, "]"))
    loginfo(str_c("Normal: [", Overweight_Obese, " <= BMI]"))
  }
  loginfo(str_c("Generate columns relating to influenza: ", Return_Influenza))
  loginfo(str_c("Generate columns relating to smoking history: ", Return_Smoker))
  loginfo(str_c("Generate columns realting to blood pressure: ", Return_Blood_Pressure))
  
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_phy_abn))
  path_phy_abn <- str_c(path_phy_abn, "Phy_")
  loginfo(str_c("Format of abnormality files: ", path_phy_abn, "{abnormality_description}", output_file_ending))
  
  #########################
  # Actual Function Start #
  #########################
  loginfo("Reading file...")
  Phy <- fread(Phy_file_name)
  loginfo("Read complete")
  rm(Phy_file_name)
  Phy <- Phy %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  
  DF_to_fill <- HealthHistoryHelper_CoreFunctionality(Physical_DF = Phy,
                                                      Main_DF = DF_to_fill,
                                                      Path_Abnormality = path_phy_abn,
                                                      File_Suffix = output_file_ending,
                                                      Gen_BMI = Return_BMI,
                                                      Und_Norm = Underweight_Normal,
                                                      Norm_Over = Normal_Overweight,
                                                      Over_Ob = Overweight_Obese,
                                                      Gen_Flu = Return_Influenza,
                                                      Gen_Smoke = Return_Smoker,
                                                      Gen_BP = Return_Blood_Pressure)
  loginfo("Processing complete")
  return(DF_to_fill)
}

#' Read in file with diagnosis information and add relevant columns to existing data frame
#' Information is limited to all before or all after
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_phy_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Phy_abnormalities/)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param Return_BMI Do you want information relating to BMI? (default = TRUE)
#' @param Underweight_Normal The value that splits BMI categorization underweight and normal (default = config$BMI_params$Underweight_Normal)
#' @param Normal_Overweight The value that splits BMI categorization normal and overweight (default = config$BMI_params$Normal_Overweight)
#' @param Overweight_Obese The value that splits BMI categorization overweight and obese (default = config$BMI_params$Overweight_Obese)
#' @param Return_Influenza Do you want information relating to influenza? (default = TRUE)
#' @param Return_Smoker Do you want information relating to smoking history? (default = TRUE)
#' @param Return_Blood_Pressure Do you want information relating to blood pressure? (default = FALSE)
#' @param Date_Column The character value representing an existing column with date information
#' @param restrict_to_before_cutoff Do you want to restrict information to before [min - cutoff; TRUE] or after [cutoff - max; FALSE] (Default = TRUE)
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_physical_date(Date_Column = Plasma_Date)
#' process_physical_date(Date_Column = Plasma_Date, Return_BMI = FALSE, Return_Blood_Pressure = TRUE)
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
                                  Return_Blood_Pressure = FALSE,
                                  Date_Column,
                                  restrict_to_before_cutoff = FALSE){
  loginfo("Processing Health History & Physical Findings...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(Date_Column)) {
    logerror("No column with Date information was specified. Process stopped.")
    return(DF_to_fill)
  }
  Phy_file_name = str_c(input_file_header, "Phy", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Phy_file_name))
  
  loginfo(str_c("Variable used for date cutoff: ", Date_Column))
  loginfo(str_c("Restriction before cutoff or after cutoff: ", ifelse(restrict_to_before_cutoff, "before", "after")))
  cutoff_prefix <- ifelse(restrict_to_before_cutoff, "Phy_PreDate_", "Phy_PostDate_")
  
  loginfo("Generate columns relating to BMI: ", Return_BMI)
  if (Return_BMI){
    loginfo(str_c("Underweight: [BMI < ", Underweight_Normal, "]"))
    loginfo(str_c("Normal: [", Underweight_Normal, " <= BMI < ", Normal_Overweight, "]"))
    loginfo(str_c("Normal: [", Normal_Overweight, " <= BMI < ", Overweight_Obese, "]"))
    loginfo(str_c("Normal: [", Overweight_Obese, " <= BMI]"))
  }
  loginfo(str_c("Generate columns relating to influenza: ", Return_Influenza))
  loginfo(str_c("Generate columns relating to smoking history: ", Return_Smoker))
  loginfo(str_c("Generate columns realting to blood pressure: ", Return_Blood_Pressure))
  
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_phy_abn))
  path_phy_abn <- str_c(path_phy_abn, cutoff_prefix)
  loginfo(str_c("Format of abnormality files: ", path_phy_abn, "{abnormality_description}", output_file_ending))
  
  #########################
  # Actual Function Start #
  #########################
  loginfo("Reading file...")
  Phy <- fread(Phy_file_name)
  loginfo("Read complete")
  rm(Phy_file_name)
  Phy <- Phy %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  # Restrict ids to only the ones in the reduced id list
  Phy <- Phy %>% filter(EMPI %in% DF_to_fill$EMPI)
  # If date variable includes time, cut off othe time
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
  # Restrict to only before or only after
  Phy <- left_join(Phy, EMPI_Date_Limit, by = "EMPI")
  rm(EMPI_Date_Limit)
  time_str = ifelse(restrict_to_before_cutoff, "before", "after")
  loginfo(str_c("Restricting Phy to ", time_str, " cutoff dates..."))
  
  Original_row_length <- nrow(Phy)
  if (restrict_to_before_cutoff){
    Phy <- Phy %>% filter(Date <= Cutoff_Date)
  } else {
    Phy <- Phy %>% filter(Date >= Cutoff_Date)
  }
  New_row_length <- nrow(Phy)
  logdebug(str_c(Original_row_length, " Phy reduced to ", New_row_length, " Phy"))
  rm(time_str, Original_row_length, New_row_length)
  DF_to_fill <- HealthHistoryHelper_CoreFunctionality(Physical_DF = Phy,
                                                      Main_DF = DF_to_fill,
                                                      Path_Abnormality = path_phy_abn,
                                                      File_Suffix = output_file_ending,
                                                      Gen_BMI = Return_BMI,
                                                      Und_Norm = Underweight_Normal,
                                                      Norm_Over = Normal_Overweight,
                                                      Over_Ob = Overweight_Obese,
                                                      Gen_Flu = Return_Influenza,
                                                      Gen_Smoke = Return_Smoker,
                                                      Gen_BP = Return_Blood_Pressure)
  loginfo("Processing complete")
  return(DF_to_fill)
}

#' Read in file with diagnosis information and add relevant columns to existing data frame
#' Information is limited to a window [min - max]
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_phy_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Phy_abnormalities/)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param Return_BMI Do you want information relating to BMI? (default = TRUE)
#' @param Underweight_Normal The value that splits BMI categorization underweight and normal (default = config$BMI_params$Underweight_Normal)
#' @param Normal_Overweight The value that splits BMI categorization normal and overweight (default = config$BMI_params$Normal_Overweight)
#' @param Overweight_Obese The value that splits BMI categorization overweight and obese (default = config$BMI_params$Overweight_Obese)
#' @param Return_Influenza Do you want information relating to influenza? (default = TRUE)
#' @param Return_Smoker Do you want information relating to smoking history? (default = TRUE)
#' @param Return_Blood_Pressure Do you want information relating to blood pressure? (default = FALSE)
#' @param min_dates The character value representing an existing column with date information
#' @param max_dates The character value representing an existing column with date information
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_physical_set_range(min_dates = "Plasma_Date_First", max_dates = "Plasma_Date_Most_Recent")
#' process_physical_set_range(min_dates = "Plasma_Date_First", max_dates = "Plasma_Date_Most_Recent", Return_BMI = FALSE, Return_Blood_Pressure = TRUE)
process_physical_set_range <- function(DF_to_fill = All_merged,
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
                                       Return_Blood_Pressure = FALSE,
                                       min_dates,
                                       max_dates){
  loginfo("Processing Health History & Physical Findings...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(min_dates)) {
    logerror("No lower bound date variable provided. Process stopped.")
    return(DF_to_fill)
  }
  if (missing(max_dates)) {
    logerror("No upper bound date variable provided. Process stopped.")
    return(DF_to_fill)
    Phy_file_name = str_c(input_file_header, "Phy", input_file_ending)
  }
  loginfo(str_c("Name of file to read in: ", Phy_file_name))
  
  loginfo("Generate columns relating to BMI: ", Return_BMI)
  if (Return_BMI){
    loginfo(str_c("Underweight: [BMI < ", Underweight_Normal, "]"))
    loginfo(str_c("Normal: [", Underweight_Normal, " <= BMI < ", Normal_Overweight, "]"))
    loginfo(str_c("Normal: [", Normal_Overweight, " <= BMI < ", Overweight_Obese, "]"))
    loginfo(str_c("Normal: [", Overweight_Obese, " <= BMI]"))
  }
  loginfo(str_c("Generate columns relating to influenza: ", Return_Influenza))
  loginfo(str_c("Generate columns relating to smoking history: ", Return_Smoker))
  loginfo(str_c("Generate columns realting to blood pressure: ", Return_Blood_Pressure))
  
  if (!dir.exists(path_phy_abn)) {dir.create(path_phy_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_phy_abn))
  path_phy_abn <- str_c(path_phy_abn, "Phy_", min_dates, "_to_", max_dates, "_")
  loginfo(str_c("Format of abnormality files: ", path_phy_abn, "{abnormality_description}", output_file_ending))
  
  #########################
  # Actual Function Start #
  #########################
  loginfo("Reading file...")
  Phy <- fread(Phy_file_name)
  loginfo("Read complete")
  rm(Phy_file_name)
  Phy <- Phy %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  
  loginfo(str_c("Restricting Phy to specified date range... "))
  Date_Range_DF <- DF_to_fill %>% select(EMPI, min_dates, max_dates)
  Phy <- right_join(Phy, Date_Range_DF, by = "EMPI")
  Original_row_length <- nrow(Phy)
  Phy <- Phy %>% filter(Date >= get(min_dates)) %>% filter(Date <= get(max_dates))
  New_row_length <- nrow(Phy)
  logdebug(str_c(Original_row_length, " Phy reduced to ", New_row_length, " Phy"))
  rm(Date_Range_DF, Original_row_length, New_row_length)
  
  DF_to_fill <- HealthHistoryHelper_CoreFunctionality(Physical_DF = Phy,
                                                      Main_DF = DF_to_fill,
                                                      Path_Abnormality = path_phy_abn,
                                                      File_Suffix = output_file_ending,
                                                      Gen_BMI = Return_BMI,
                                                      Und_Norm = Underweight_Normal,
                                                      Norm_Over = Normal_Overweight,
                                                      Over_Ob = Overweight_Obese,
                                                      Gen_Flu = Return_Influenza,
                                                      Gen_Smoke = Return_Smoker,
                                                      Gen_BP = Return_Blood_Pressure)
  loginfo("Processing complete")
  return(DF_to_fill)
}
