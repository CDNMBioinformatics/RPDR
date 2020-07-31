require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)
require(readr)
require(zeallot) # %<-% (multiple variable assignment)
require(sqldf)

####################################################################################################
######################################## General functions  ########################################
####################################################################################################
# If query came from Biobank, use this as step one
process1_biobankIDs <- function(input_file_header = config$rpdr_file_header,
                                input_file_ending = config$rpdr_file_ending){
  loginfo("Processing biobank ids file... ")
  BiobankIDs <- data.table(fread(str_c(input_file_header, "Bib", input_file_ending)))
  BiobankIDs <- BiobankIDs %>% select(Subject_Id, EMPI) %>% rename(Biobank_Subject_ID = Subject_Id)
  loginfo(str_c(nrow(BiobankIDs), " subjects processed"))
  return(BiobankIDs)
}

# If query came from RPDR, use this as step one
process1_demographics <- function(input_file_header = config$rpdr_file_header,
                                  input_file_ending = config$rpdr_file_ending){
  loginfo("Processing demographics file...")
  Demographics <- data.table(fread(str_c(input_file_header, "Dem", input_file_ending)))
  Demographics <- Demographics %>% select(EMPI, Gender, Race, Date_of_Birth, Age, Date_Of_Death, Vital_status)
  loginfo(str_c(nrow(Demographics), " subjects processed"))
  return(Demographics)
}

process2_demographics <- function(DF_to_fill = All_merged,
                                  input_file_header = config$rpdr_file_header,
                                  input_file_ending = config$rpdr_file_ending){
  loginfo("Processing demographics file...")
  Demographics <- data.table(fread(str_c(input_file_header, "Dem", input_file_ending)))
  Demographics <- Demographics %>% select(EMPI, Gender, Race, Date_of_Birth, Age, Date_Of_Death, Vital_status)
  DF_to_fill <- left_join(DF_to_fill, Demographics, by = "EMPI")
  loginfo("Gender, race,  date of birth, age, date of death, and vital status information have been added")
  rm(Demographics)
  return(DF_to_fill)
}

process2_biobankIDs <- function(DF_to_fill = All_merged,
                                input_file_header = config$rpdr_file_header,
                                input_file_ending = config$rpdr_file_ending){
  loginfo("Processing biobank ids file... ")
  BiobankIDs <- data.table(fread(str_c(input_file_header, "Bib", input_file_ending)))
  BiobankIDs <- BiobankIDs %>% select(Subject_Id, EMPI) %>% rename(Biobank_Subject_ID = Subject_Id)
  DF_to_fill <- left_join(DF_to_fill, BiobankIDs, by = "EMPI")
  loginfo("Available BiobankIds have been added")
  rm(BiobankIDs)
  return(DF_to_fill)
}

process2_identifable_information <- function(DF_to_fill = All_merged,
                                   input_file_header = config$rpdr_file_header,
                                   input_file_ending = config$rpdr_file_ending){
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
  DF_to_fill <- left_join(DF_to_fill, Identifiable, by = "EMPI")
  rm(Identifiable)
  return(DF_to_fill)
}

####################################################################################################
#################################### Deidentification functions ####################################
####################################################################################################
process3_deidentified <- function(DF_to_fill = All_merged,
                                  input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  Deidentified <- Deidentified %>% select(`Biobank Subject ID`, contains("Asthma"))
  names(Deidentified) <- gsub("(.*)_current_.*", "\\1",
                              gsub("(.*)_no_history_.*", "Non\\1",
                                   gsub("_+", "_",
                                        gsub("( |\\(|\\)|-|/|\\[|\\])", "_",
                                             gsub("(.* )(\\[.*)", "\\1",
                                                  names(Deidentified))))))
  Deidentified <- Deidentified %>%
    mutate_at(vars(contains("Asthma")), ~ifelse(.x == "Yes", 1, 0)) %>%
    mutate(Biobank_Subject_ID = as.integer(Biobank_Subject_ID))
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and NonAsthma identifiers have been added")
  rm(Deidentified)
  return(DF_to_fill)
}

process3_deidentified_asthma_copd <- function(DF_to_fill = All_merged,
                                              input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  Deidentified <- Deidentified %>% 
    select(`Biobank Subject ID`, contains("Asthma"), contains("COPD"))
  names(Deidentified) <- gsub("(.*)_current_.*", "\\1",
                              gsub("(.*)_no_history_.*", "Non\\1",
                                   gsub("_+", "_",
                                        gsub("( |\\(|\\)|-|/|\\[|\\])", "_",
                                             gsub("(.* )(\\[.*)", "\\1",
                                                  names(Deidentified))))))
  Deidentified <- Deidentified %>%
    mutate_at(vars(contains("Asthma"), contains("COPD")), ~ifelse(.x == "Yes", 1, 0)) %>%
    select(Biobank_Subject_ID, Asthma, NonAsthma, COPD, NonCOPD) %>%
    mutate(Biobank_Subject_ID = as.integer(Biobank_Subject_ID))
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and COPD identifiers have been added")
  rm(Deidentified)
  return(DF_to_fill)
}

process3_deidentified_asthma_tobacco <- function(DF_to_fill = All_merged,
                                                 input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  Deidentified <- Deidentified %>% 
    select(`Biobank Subject ID`, contains("Asthma"), contains("Tobacco"))
  names(Deidentified) <- gsub("Tobacco_User_Never", "NonSmoker", 
                              gsub("_$", "",
                                   gsub("(.*)_current_.*", "\\1",
                                        gsub("(.*)_no_history_.*", "Non\\1",
                                             gsub("_+", "_",
                                                  gsub("( |\\(|\\)|-|/|\\[|\\])", "_",
                                                       gsub("(.* )(\\[.*)", "\\1",
                                                            names(Deidentified))))))))
  Deidentified <- Deidentified %>%
    mutate_at(vars(contains("Asthma"), contains("NonSmoker")), ~ifelse(.x == "Yes", 1, 0)) %>%
    select(Biobank_Subject_ID, Asthma, NonAsthma, everything()) %>%
    mutate(Biobank_Subject_ID = as.integer(Biobank_Subject_ID))
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and COPD identifiers have been added")
  rm(Deidentified)
  return(DF_to_fill)
}

process_List <- function(List_df = Deidentified){
  loginfo("Processing List All Values items...")
  List_df <- List_df %>% mutate_at(vars(contains("List")), ~gsub("\\[((\\d|\\.)+)\\]", "\\1", .x))
  List_df <- List_df %>% mutate_at(vars(contains("List")), ~gsub(" ", ";", .x))
  List_df <- List_df %>% mutate_at(vars(contains("List")), ~gsub(";$", "", .x))
  List_df <- List_df %>% mutate_at(vars(contains("List")), ~gsub("^$", NA, .x))
  return(List_df)
}

process_dates <- function(DF = Deidentified){
  loginfo("Processing date ranges...")
  DF <- DF %>% mutate_at(vars(contains("Date")), ymd_hms)
  Date_pairs <- gsub("^(.*_)Date_First$", "\\1", DF %>% select(contains("Date_First")) %>% names())
  for (Group in Date_pairs){
    DF <- DF %>% mutate(!!(as.symbol(str_c(Group, "Range_Years"))) := 
                          time_length(interval(!!(as.symbol(str_c(Group, "Date_First"))),
                                               !!(as.symbol(str_c(Group, "Date_Most_Recent")))),
                                      unit = "year")) %>% 
      select(Biobank_Subject_ID:matches(str_c("^", Group, "Date_Most_Recent")),
             str_c(Group, "Range_Years"), everything())
  }
  return(DF)
}
process3_deidentified_complex <- function(DF_to_fill = All_merged,
                                          input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  names(Deidentified) <- gsub("^(.*)_$", "\\1",
                              gsub("_+", "_",
                                   gsub(" |-|\\(|\\)|\\[|\\]|;|:|/|&", "_", names(Deidentified))))
  Deidentified <- Deidentified %>%
    rename(Asthma = Asthma_current_or_past_history_PPV_0.90_Existence_Yes_No,
           NonAsthma = Asthma_no_history_NPV_0.99_Existence_Yes_No) %>%
    mutate(Asthma = ifelse(Asthma == "Yes", 1, 0), NonAsthma = ifelse(NonAsthma == "Yes", 1, 0))
  Deidentified <- Deidentified %>% select(-(Gender:Vital_Status)) %>%
    select(-(starts_with("BMI"))) %>%
    mutate(Biobank_Subject_ID = as.integer(Biobank_Subject_ID))
  Deidentified <- process_List(Deidentified)
  Deidentified <- process_dates(Deidentified)
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and NonAsthma, BMI, and Spirometry identifiers have been added")
  rm(Deidentified)
  return(DF_to_fill)
}
process3_deidentified_complex2 <- function(DF_to_fill = All_merged,
                                           input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  names(Deidentified) <- gsub("^(.*)_$", "\\1",
                              gsub("_+", "_",
                                   gsub(" |-|\\(|\\)|\\[|\\]|;|:|/|&", "_", names(Deidentified))))
  Deidentified <- Deidentified %>%
    rename(Asthma_0.95PPV = `Asthma_current_or_past_history_custom_PPV_>=_0.95PPV_Existence_Yes_No`,
           COPD_0.94PPV = `COPD_current_or_past_history_custom_PPV_>=_0.94PPV_Existence_Yes_No`,
           T2DM_0.99PPV = T2DM_current_or_past_history_PPV_0.99_Existence_Yes_No,
           Obesity_0.97PPV = `Obesity_current_or_past_history_custom_PPV_>=_0.97PPV_Existence_Yes_No`) %>%
    mutate_at(vars(contains("PPV")), ~ifelse(.x == "Yes", 1, 0))
  Deidentified <- process_List(Deidentified)
  Deidentified <- Deidentified %>% select(Biobank_Subject_ID, Body_Mass_Index_BMI_Existence_Yes_No:Obesity_0.97PPV)
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  rm(Deidentified)
  return(DF_to_fill)
}

process3_deidentified_ppv_only <- function(DF_to_fill = All_merged,
                                           input_file_name = config$biobank_file_name){
  loginfo("Processing biobank file...")
  Deidentified <- fread(input_file_name)
  names(Deidentified) <- gsub("^(.*)_$", "\\1",
                              gsub("_+", "_",
                                   gsub(" |-|\\(|\\)|\\[|\\]|;|:|/|&", "_", names(Deidentified))))
  Deidentified <- Deidentified %>% select(Biobank_Subject_ID, contains("PPV"))
  names(Deidentified) <- gsub("^(.*)_current_or_past_history_custom_PPV_>=_((\\d|\\.)*PPV)_(.*)$",
                              "\\1_\\2_\\4",
                              names(Deidentified))
  Deidentified <- Deidentified %>% mutate_at(vars(contains("List_of_All_Values")),
                                             ~as.numeric(gsub("\\[(.*)\\]", "\\1", .x)),
                                             vars(contains("Existence")),
                                             ~ifelse(grepl("Yes", .x), 1, 0))
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  rm(Deidentified)
  return(DF_to_fill)
}

####################################################################################################
####################################### Diagnosis functions  #######################################
####################################################################################################
process4_diagnoses <- function(DF_to_fill = All_merged,
                               input_file_header = config$rpdr_file_header,
                               input_file_ending = config$rpdr_file_ending,
                               path_dia_abn = str_c(config$data_dir, "Diagnoses_abnormalities/"),
                               Diagnoses_Of_Interest,
                               Group_Info = TRUE,
                               Individual_Info = TRUE,
                               write_files = config$create_intermediates,
                               output_file_header = config$intermediate_files_dir,
                               output_file_ending = config$general_file_ending){
  if (missing(Diagnoses_Of_Interest)){
    logerror("No list of Diagnoses were specified. Process stopped.")
    return(DF_to_fill)
  }
  loginfo("Processing diagnoses file...")
  Diagnoses <- data.table(fread(str_c(input_file_header, "Dia", input_file_ending))) %>% 
    arrange(EMPI, Date)
  if (!dir.exists(path_dia_abn)) {dir.create(path_dia_abn)}
  # Get the "Any exist" first to lower the search group/increase speed later
  for (Grouping_Name in names(Diagnoses_Of_Interest)){
    Group_Header = str_c("Any_", gsub(" ", "_", Grouping_Name), "_diagnosis")
    Group <- Diagnoses %>%
      filter(grepl(str_c(Diagnoses_Of_Interest[[Grouping_Name]], collapse = "|"), Diagnosis_Name))
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Diagnosis_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects have any ", Grouping_Name, " diagnosis"))
    rm(Output_Columns)
    
    if (Individual_Info){
      # Look for the individual diagnoses
      for (Diagnosis in Diagnoses_Of_Interest[[Grouping_Name]]){
        Subgroup_Header <- gsub(" ", "_", Diagnosis)
        Subgroup <- Group %>% filter(Diagnosis_Name == Diagnosis) %>% group_by(EMPI)
        Dia_abn <- Subgroup[duplicated(Subgroup) | duplicated(Subgroup, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " completely duplicated row(s) out of ",
                        nrow(Subgroup), " found. Duplicates removed."))
          fwrite(Dia_abn, str_c(path_dia_abn, "Abnormality_1_Duplicate_rows_",
                                Subgroup_Header, output_file_ending))
          Subgroup <- Subgroup %>% unique()
        }
        Subgroup2 <- Subgroup %>% select(EMPI, Date, Diagnosis_Name, Hospital)
        Dia_abn <- Subgroup[duplicated(Subgroup2) | duplicated(Subgroup2, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Subgroup),
                        " found. ", sum(duplicated(Subgroup2, fromLast = TRUE)), " rows removed."))
          fwrite(Dia_abn, str_c(path_dia_abn, "Abnormality_2_Duplicate_rows_",
                                Subgroup_Header, output_file_ending))
          Subgroup <- distinct(Subgroup, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup2)
        
        Output_Columns <- Subgroup %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Date, collapse = ";"))
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Diagnosis), " diagnosis"))
        
        if(write_files){
          fwrite(Subgroup, str_c(output_file_header, Subgroup_Header, output_file_ending))
        }
        rm(Subgroup_Header)
        rm(Subgroup, Output_Columns, Dia_abn)
      }
    }
    Group <- Group %>% unique()
    Group2 <- Group %>% select(EMPI, Date, Diagnosis_Name, Hospital)
    Dia_abn <- Group[duplicated(Group2) | duplicated(Group2, fromLast=TRUE),]
    if (nrow(Dia_abn) > 0){
      logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Group),
                    " found. ", sum(duplicated(Group2, fromLast = TRUE)), " rows removed."))
      fwrite(Dia_abn, str_c(path_dia_abn, "Abnormality_2_Duplicate_rows_",
                            Group_Header, output_file_ending))
      Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
    }
    rm(Group2)
    Output_Columns <- Group %>% group_by(EMPI) %>% mutate(Date = mdy(Date)) %>% arrange(Date) %>%
      summarise(!!(as.symbol(str_c(Group_Header, "_total_dates"))) := n(),
                !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    
    if(write_files){
      fwrite(Group, str_c(output_file_header, Group_Header, output_file_ending))
    }
    rm(Group_Header)
    rm(Group, Output_Columns, Dia_abn)
  }
  
  rm(Diagnoses)
  return(DF_to_fill)
}


####################################################################################################
####################################### Medication functions #######################################
####################################################################################################

process5_medications <- function(DF_to_fill = All_merged,
                                 input_file_header = config$rpdr_file_header,
                                 input_file_ending = config$rpdr_file_ending,
                                 path_med_abn = str_c(config$data_dir, "Medication_abnormalities/"),
                                 Medications_Of_Interest,
                                 Group_Column = config$medication_group,
                                 Individual_Info = TRUE,
                                 Group_Info = TRUE,
                                 Daily_Dose_Info = FALSE,
                                 write_files = config$create_intermediates,
                                 output_file_header = config$intermediate_files_dir,
                                 output_file_ending = config$general_file_ending,
                                 merged_group_name,
                                 nebs = FALSE){
  if (missing(Medications_Of_Interest)){
    logerror("No list of Medications were specified. Process stopped.")
    return(DF_to_fill)
  }
  loginfo("Processing medications file...")
  Medications <- fread(str_c(input_file_header, "Med", input_file_ending))
  if (!dir.exists(path_med_abn)) {dir.create(path_med_abn)}
  
  loginfo("Applying Med_Map...")
  Med_Map_df <- fread(str_c(general_path, "/Medication_Mapping.txt"))
  logdebug("Available folders and number of medication groups:")
  logdebug(t(Med_Map_df %>% group_by(Medication_Biobank_Folder) %>% summarise(Count = n())))
  # Select relevant files
  Relevant_Meds <- data.frame(Medication_Name = character(),
                              Medication_Biobank_Folder = character(),
                              Medication_GWAS_Group = character(),
                              Medication_Pegasus_Group = character(),
                              Search_Term = character(),
                              Ignore.case = logical(),
                              Perl = logical())
  for (Grouping_Name in names(Medications_Of_Interest)){
    if (length(Medications_Of_Interest[[Grouping_Name]]) == 0){
      Mini_Med_Map <- Med_Map_df %>% filter(!!(as.symbol(Group_Column)) %in% Grouping_Name)
      Medications_Of_Interest[[Grouping_Name]] = Mini_Med_Map %>% pull(Medication_Name)
    } else {
      Mini_Med_Map <- Med_Map_df %>% filter(!!(as.symbol(Group_Column)) %in% Grouping_Name,
                                            Medication_Name %in% Medications_Of_Interest[[Grouping_Name]])
    }
    if (nrow(Mini_Med_Map) == 0){
      logerror("Med grouping does not exist. Please try again.")
      return(DF_to_fill)
    }
    Relevant_Meds <- rbind(Relevant_Meds, Mini_Med_Map)
  }
  rm(Mini_Med_Map, Med_Map_df)
  Med_Row_ids <- 1:nrow(Relevant_Meds)
  Medications_condensed <- Medications %>% group_by(Medication) %>% summarise()
  Med_all_subgroups <- data.frame(Medication = character(),
                                  Medication_Name = character(),
                                  Medication_Biobank_Folder = character(),
                                  Medication_GWAS_Group = character(),
                                  Medication_Pegasus_Group = character())
  for (mri in Med_Row_ids){
    Med_subgroup <- Medications_condensed %>% filter(grepl(Relevant_Meds[mri,]$Search_Term,
                                                           Medication,
                                                           ignore.case = Relevant_Meds[mri,]$Ignore.case,
                                                           perl = Relevant_Meds[mri,]$Perl)) %>%
      mutate(Medication_Name  = Relevant_Meds[mri,]$Medication_Name,
             Medication_Biobank_Folder = Relevant_Meds[mri,]$Medication_Biobank_Folder,
             Medication_GWAS_Group = Relevant_Meds[mri,]$Medication_GWAS_Group,
             Medication_Pegasus_Group = Relevant_Meds[mri,]$Medication_Pegasus_Group)
    Med_all_subgroups <- rbind(Med_all_subgroups, Med_subgroup)
    rm(Med_subgroup)
  }
  Medications <- right_join(Medications, Med_all_subgroups, by = "Medication")
  rm(Med_all_subgroups, mri, Medications_condensed, Relevant_Meds, Med_Row_ids)
  # Med Cleanup
  Medications <- Medications %>% mutate(Medication_Date = mdy(Medication_Date)) %>% arrange(EMPI, Medication_Date)
  # Clean up Additional_Info and generate Daily Dosage and Notes
  # - Get MCG units (if ML: get MG from medication; if MG: convert to MCG) or PUFF count
  # - Get FREQ
  # - Note PRN
  # - Calculate Daily Dosages
  # - write Notes
  Medications <- Medications %>% mutate(Additional_Info = gsub("PUFFS", "PUFF", Additional_Info, ignore.case = TRUE),
                                        Additional_Info = gsub("mcg", "MCG", Additional_Info, ignore.case = TRUE),
                                        Additional_Info = gsub("ml", "ML", Additional_Info, ignore.case = TRUE),
                                        Additional_Info = gsub(" of fluti", "", Additional_Info),
                                        Additional_Info = gsub("Inhl", "INH", Additional_Info, ignore.case = TRUE),
                                        Additional_Info = gsub("Nebu", "NEB", Additional_Info, ignore.case = TRUE),
                                        DOSE_MCG = as.numeric(gsub("^.*DOSE=((\\d|\\.)+) MCG.*$", "\\1", Additional_Info)),
                                        DOSE_MG = as.numeric(gsub("^.*DOSE=((\\d|\\.)+) MG.*$", "\\1", Additional_Info)),
                                        DOSE_ML = as.numeric(gsub("^.*DOSE=((\\d|\\.)+) ML.*$", "\\1", Additional_Info)),
                                        DOSE_MG_ML = ifelse(is.na(DOSE_ML), NA, Medication),
                                        DOSE_MG_ML_mg = as.numeric(gsub("^.* ((\\d|\\.)+) *mg.*$", "\\1", DOSE_MG_ML, ignore.case = TRUE)),
                                        DOSE_MG_ML_ml = as.numeric(gsub("^.*/ *((\\d|\\.)*) *ml.*$", "\\1", DOSE_MG_ML, ignore.case = TRUE)),
                                        DOSE_MG_ML_ml = ifelse(!is.na(DOSE_ML) & is.na(DOSE_MG_ML_ml), 1, DOSE_MG_ML_ml),
                                        DOSE_MG = ifelse(!(is.na(DOSE_ML)), DOSE_ML*DOSE_MG_ML_mg/DOSE_MG_ML_ml, DOSE_MG),
                                        DOSE_MCG = ifelse(!(is.na(DOSE_MG)), DOSE_MG*1000, DOSE_MCG),
                                        DOSE_PUFF = as.numeric(gsub("^.*DOSE=((\\d|\\.)*) (INHALATION|PUFF|SPRAY).*$", "\\1", Additional_Info)),
                                        FREQ = NA,
                                        FREQ = ifelse(grepl("Daily|Nightly|Once|QAM|QD|QNOON|QPM|Q24", Additional_Info), 1, FREQ),
                                        FREQ = ifelse(grepl("BID|Q12H", Additional_Info), 2, FREQ),
                                        FREQ = ifelse(grepl("QAC|TID|Q8H", Additional_Info), 3, FREQ),
                                        FREQ = ifelse(grepl("4x Daily|QID|Q6H", Additional_Info), 4, FREQ),
                                        FREQ = ifelse(grepl("Q4H", Additional_Info), 6, FREQ),
                                        FREQ = ifelse(grepl("Q3H", Additional_Info), 8, FREQ),
                                        FREQ = ifelse(grepl("Q2H", Additional_Info), 12, FREQ),
                                        FREQ = ifelse(grepl("Q1H", Additional_Info), 24, FREQ),
                                        PRN = ifelse(grepl("PRN", Additional_Info), TRUE, FALSE),
                                        DAILY_DOSE_MCG = ifelse(is.na(DOSE_MCG), NA, DOSE_MCG * ifelse(is.na(FREQ), 1, FREQ)),
                                        DAILY_DOSE_PUFF = ifelse(is.na(DOSE_PUFF), NA, DOSE_PUFF * ifelse(is.na(FREQ), 1, FREQ)),
                                        DAILY_DOSE = ifelse(!is.na(DAILY_DOSE_MCG), str_c("MCG: ", DAILY_DOSE_MCG), NA),
                                        DAILY_DOSE = ifelse(!is.na(DAILY_DOSE_PUFF), str_c("Puffs: ", DAILY_DOSE_PUFF), DAILY_DOSE),
                                        DAILY_DOSE = ifelse(is.na(DAILY_DOSE) & !is.na(FREQ), str_c("FREQ: ", FREQ), DAILY_DOSE),
                                        NOTES = ifelse(PRN, "Prescribed as needed", ""),
                                        NOTES = ifelse(is.na(DAILY_DOSE) | grepl("FREQ", DAILY_DOSE), str_c(NOTES, ifelse(NOTES == "", "", "; "), "Unknown dosage"), NOTES), 
                                        NOTES = ifelse(is.na(FREQ), str_c(NOTES, ifelse(NOTES == "", "", "; "), "Unknown frequency"), NOTES)) %>%
    select(-c(DOSE_MCG, DOSE_MG, DOSE_ML, DOSE_MG_ML, DOSE_MG_ML_mg, DOSE_MG_ML_ml, DOSE_PUFF, FREQ, PRN))
  
  if ("Medication_Date_Detail" %in% names(Medications)){
    Med_abn <- Medications %>% filter(Medication_Date_Detail == "Removed")
    if (nrow(Med_abn) > 0){
      logwarn(str_c(nrow(Med_abn), " row(s) out of ", nrow(Medications), " have been removed due having the flag 'Removed'."))
      fwrite(Med_abn, str_c(path_med_abn, "Abnormality_1_Removed_Flag_All_Medications", output_file_ending))
      Medications <- Medications %>% filter(Medication_Date_Detail != "Removed")
    }
  }
  loginfo("Creating medication output columns...")
  # Get the "Any exist" first to lower the search group/increase speed later
  if (!missing(merged_group_name)){
    Merged_Group_Header = str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", merged_group_name)))
    Output_Columns = Medications %>% group_by(EMPI) %>% select(EMPI, Medication_Name) %>%
      unique() %>% summarise(!!(as.symbol(Merged_Group_Header)) := "Yes")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(!!(as.symbol(Merged_Group_Header)) := ifelse(is.na(!!(as.symbol(Merged_Group_Header))),
                                                          "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", merged_group_name))
  }
  for (Grouping_Name in names(Medications_Of_Interest)){
    Group_Header = str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)))
    Group <- Medications %>%
      filter(grepl(Grouping_Name, !!(as.symbol(Group_Column))),
             Medication_Name %in% Medications_Of_Interest[[Grouping_Name]])
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Medication_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes")
    if (Group_Info){
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      DF_to_fill <- DF_to_fill %>%
        mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
    }
    loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", Grouping_Name))
    rm(Output_Columns)
    if(Individual_Info){
      # Look for the individual prescriptions
      for (Med_Name in Medications_Of_Interest[[Grouping_Name]]){
        Subgroup_Header <- str_c(gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)),
                                 "_", gsub("/| ", "_", Med_Name))
        Subgroup <- Group %>% filter(Medication_Name == Med_Name) %>% group_by(EMPI)
        Med_abn <- Subgroup[duplicated(Subgroup) | duplicated(Subgroup, fromLast = TRUE),]
        if (nrow(Med_abn) > 0){
          logwarn(str_c(nrow(Med_abn), " completely duplicated row(s) out of ", nrow(Subgroup), " found. Duplicates removed."))
          fwrite(Med_abn, str_c(path_med_abn, "Abnormality_2_Duplicate_rows_", Subgroup_Header, output_file_ending))
          Subgroup <- Subgroup %>% unique()
        }
        Subgroup2 <- Subgroup %>% select(EMPI, Medication_Date, Medication, Hospital)
        Med_abn <- Subgroup[duplicated(Subgroup2) | duplicated(Subgroup2, fromLast = TRUE),]
        if (nrow(Med_abn) > 0){
          logwarn(str_c(nrow(Med_abn), " partially duplicated row(s) out of ", nrow(Subgroup), " found. ",
                        sum(duplicated(Subgroup2, fromLast = TRUE)), " rows removed."))
          fwrite(Med_abn, str_c(path_med_abn, "Abnormality_3_Duplicate_rows_", Subgroup_Header, output_file_ending))
          Subgroup <- distinct(Subgroup, EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup2)
        # Start writing output
        Output_Columns <- Subgroup %>% group_by(EMPI) %>% summarise(!!(as.symbol(Subgroup_Header)) := "Yes")
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        DF_to_fill <- DF_to_fill %>%
          mutate(!!(as.symbol(Subgroup_Header)) := ifelse(is.na(!!(as.symbol(Subgroup_Header))), "No", "Yes"))
        Output_Columns <- Subgroup %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_first_date"))) := first(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_last_date"))) := last(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Medication_Date, collapse = ";"))
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        Output_Columns <- Subgroup %>% group_by(EMPI, Medication) %>% summarise(Medication_Occurances = n()) %>%
          group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                    !!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"))
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        if (Daily_Dose_Info){
          Output_Columns <- Subgroup %>% group_by(EMPI) %>%
            summarise(!!(as.symbol(str_c(Subgroup_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_notes"))) := paste(NOTES, collapse= "|"))
          DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        }
        loginfo(str_c(nrow(Output_Columns), " subjects were prescribed ", tolower(Med_Name), "."))
        logdebug(str_c("nlines ", Subgroup_Header, " after completion: ", nrow(Subgroup)))
        if(write_files){
          fwrite(Subgroup, str_c(output_file_header, Group_Header, "_", Subgroup_Header, output_file_ending))
        }
        rm(Subgroup_Header)
        rm(Subgroup, Output_Columns, Med_abn)
      }
      rm(Med_Name)
    }
    if(Group_Info){
      # Now that all the errors have been noted.. add more count/date columns
      Group <- Group %>% distinct(EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
      Output_Columns <- Group %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_total_dates"))) := n(),
                  !!(as.symbol(str_c(Group_Header, "_first_date"))) := first(Medication_Date),
                  !!(as.symbol(str_c(Group_Header, "_last_date"))) := last(Medication_Date),
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"))
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication_Name) %>% summarise(Medication_Occurances = n()) %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"))
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication) %>% summarise(Medication_Occurances = n()) %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"))
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      if (Daily_Dose_Info){
        Output_Columns <- Group %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_notes"))) := paste(NOTES, collapse= "|"))
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      }
      # Rearrange so "Any" columns are together for easier viewing/understanding
      DF_to_fill <- DF_to_fill %>% select(EMPI:Group_Header, starts_with(Group_Header), everything())
      rm(Output_Columns)
      if(write_files){
        fwrite(Group, str_c(output_file_header, Group_Header, output_file_ending))
      }
    }
    rm(Group_Header)
    rm(Group)
  }
  if (!missing(merged_group_name)){
    Medications_distinct <- Medications %>% distinct(EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
    Merged_Group_Header = str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", merged_group_name)))
    Output_Columns <- Medications_distinct %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_total_dates"))) := n(),
                !!(as.symbol(str_c(Merged_Group_Header, "_first_date"))) := first(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_last_date"))) := last(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication_Name) %>% summarise(Medication_Occurances = n()) %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication) %>% summarise(Medication_Occurances = n()) %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    if (Daily_Dose_Info){
      Output_Columns <- Medications_distinct %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_notes"))) := paste(NOTES, collapse= "|"))
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    }
    if (nebs){
      Output_Columns <- Medications_distinct %>% filter(grepl("[Nn]eb", Medication)) %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) := "Yes")
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      DF_to_fill <- DF_to_fill %>%
        mutate(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) :=
                 ifelse(is.na(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer")))),
                        "No", "Yes"))
    }
    # Rearrange so "Any" columns are together for easier viewing/understanding
    DF_to_fill <- DF_to_fill %>% select(EMPI:Merged_Group_Header, starts_with(Merged_Group_Header), everything())
    rm(Output_Columns)
    if(write_files){
      fwrite(Medications_distinct, str_c(output_file_header, Merged_Group_Header, output_file_ending))
    }
    rm(Medications_distinct)
  }
  return(DF_to_fill)
}

####################################################################################################
########################################## Lab functions  ##########################################
####################################################################################################
Create_ACTH_Cortisol_DHEA_Output_Columns <- function(ACTH_Cortisol_DHEA_Group,
                                                     Group_Header,
                                                     DF_to_fill){
  # Splitting up Output Columns step into substeps because some median date selection doesn't work well with summarise
  OC_pregroup <- ACTH_Cortisol_DHEA_Group %>% group_by(EMPI, Reference_Units, Seq_Date, Seq_Time) %>%
    summarise(nResults_Per_Time = n(),
              MinResult = min(Result),
              MaxResult = max(Result),
              MedianResult = median(Result),
              MeanResult = mean(Result)) %>%
    arrange(Seq_Time) %>% group_by(EMPI, Reference_Units, Seq_Date) %>%
    summarise(nResults_Per_Date = sum(nResults_Per_Time),
              nTimes_Per_Date = n(),
              FirstMin = first(MinResult),
              FirstMax = first(MaxResult),
              FirstMedian = first(MedianResult),
              FirstMean = first(MeanResult),
              FirstTime = first(Seq_Time)) %>%
    mutate(Seq_Date_Time = ifelse(is.na(FirstTime), as.character(Seq_Date), str_c(Seq_Date, " ", FirstTime))) %>% 
    group_by(EMPI, Reference_Units) %>%
    rename(!!as.symbol(str_c(Group_Header, "_Reference_Units")) := Reference_Units)
  DF_to_fill <- left_join(DF_to_fill, OC_pregroup %>% summarise(), by = "EMPI")
  OC_pregroup <- OC_pregroup %>% group_by(EMPI) %>% arrange(Seq_Date_Time)
  
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_nTotalDates")) := n(),
              !!as.symbol(str_c(Group_Header, "_nTotalDatesTimes")) := sum(nTimes_Per_Date),
              !!as.symbol(str_c(Group_Header, "_nTotalResults")) := sum(nResults_Per_Date),
              !!as.symbol(str_c(Group_Header, "_All_Seq_Date_Times")) := paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Min_Result")) := min(FirstMin),
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_First")) := Seq_Date_Time[first(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_All_Min_Results")) := paste(FirstMin, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Max_Result")) := max(FirstMax),
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_First")) := Seq_Date_Time[first(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_All_Max_Results")) := paste(FirstMax, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>% filter(abs(median(FirstMedian) - FirstMedian) == median(FirstMedian) - FirstMedian) %>%
    arrange(EMPI, FirstMedian) %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Median_Result")) := max(FirstMedian),
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_First_or_closest_below")) := Seq_Date_Time[first(which(FirstMedian == max(FirstMedian)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_Last_or_closest_below")) := Seq_Date_Time[last(which(FirstMedian == max(FirstMedian)))])
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  # All Median is done separately because the previous computation dumps rows that we don't want dumped for the paste
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_All_Median_Results")) := paste(FirstMedian, collapse = ";"),
              !!as.symbol(str_c(Group_Header, "_Overall_Mean_Result")) := mean(FirstMean),
              !!as.symbol(str_c(Group_Header, "_All_Mean_Results")) := paste(FirstMean, collapse = ";"))
  
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  loginfo(str_c(nrow(Output_Columns), " subjects have a ", Group_Header, " test performed with useful information"))
  rm(OC_pregroup, Output_Columns)
  return(DF_to_fill)
}

process6_ACTH_labs <- function(DF_to_fill = All_merged,
                               input_file_header = config$rpdr_file_header,
                               input_file_ending = config$rpdr_file_ending,
                               path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                               skip_ACTH = config$ACTH_params$skip_ACTH,
                               strict = config$ACTH_params$strict,
                               create_cortisol_group = config$ACTH_params$create_cortisol_group,
                               write_files = config$create_intermediates,
                               output_file_header = config$intermediate_files_dir,
                               output_file_ending = config$general_file_ending){
  loginfo("Processing labs file...")
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI, Seq_Date_Time)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  
  # Only care about ACTH, Cortisol, and/or DHEA(s) in this function
  Labs <- Labs %>% filter(grepl("ACTH|Cortisol|DHEA", Group_Id))
  
  # Clean data (1):
  # (A): Change "Less than x" to "< x" and "Greater than x" to "> x"; Get rid of an extra spacing
  # (B): Only include values that are digits, decimal starts, < x, > x
  Lab_abn <- Labs %>% filter(!grepl("^ *(<|>|LESS THAN|GREATER THAN)* *(\\d|\\.)", toupper(Result)))
  logwarn(str_c(nrow(Lab_abn), " entries out of ", nrow(Labs), " removed due to missingness or corrupt result entries"))
  fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_1_Missingness_or_corrupt_result", output_file_ending))
  Labs <- Labs %>%
    mutate(Result = gsub("Less than", "<", Result, ignore.case = TRUE),
           Result = gsub("Greater than", ">", Result, ignore.case = TRUE),
           Result = gsub(" ", "", Result)) %>%
    filter(grepl("^(<|>)*(\\d|\\.)", Result))
  
  # Clean data (2)
  #   Change <0 to 0; Change all <x to x-min(x)/10; Change all >x to x + 0.11111
  Labs <- Labs %>%
    mutate(Result = gsub("<0((\\.|0)*)$", "0", Result),
           Result = gsub("(.*(\\d|\\.)+).*", "\\1", Result),
           LessThanX = as.numeric(ifelse(grepl("<", Result), gsub("<((\\d|\\.)+)", "\\1", Result), NA)),
           GreaterThanX = as.numeric(ifelse(grepl(">", Result), gsub(">((\\d|\\.)+)", "\\1", Result), NA)),
           Result = as.numeric(Result),
           Result = ifelse(is.na(LessThanX), Result, LessThanX - min(LessThanX, na.rm = TRUE) / 10),
           Result = ifelse(is.na(GreaterThanX), Result, GreaterThanX + 1/9)) %>%
    select(-c(LessThanX, GreaterThanX))
  
  # Clean data (3)
  #   Remove duplicates
  Lab_abn <- Labs[(duplicated(Labs)),]
  logwarn(str_c(nrow(Lab_abn), " completely duplicated row(s) out of ", nrow(Labs), " removed"))
  fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_2_Duplicate_rows", output_file_ending))
  Labs <- Labs %>% unique()
  
  # Clean data (4)
  # (A): Change all reference units to lower case; Change mcg (micrograms) to ug (micrograms);
  #      Find units listed in result text
  # (B): Remove mismatches between units in result text and reference units
  Labs <- Labs %>%
    mutate(Reference_Units = tolower(Reference_Units),
           Reference_Units = gsub("mc", "u", Reference_Units),
           Result_Text = gsub("mcg/", "ug/", Result_Text, ignore.case = TRUE),
           Result_Text_Units = ifelse(grepl(".*([unp]g/[dm]l).*", Result_Text, ignore.case = TRUE),
                                      gsub(".*([unp]g */ *[dm]l).*", "\\1", Result_Text, ignore.case = TRUE),
                                      ""),
           Result_Text_Units = tolower(Result_Text_Units))
  Lab_abn <- Labs %>% filter(Result_Text_Units != "" & Result_Text_Units != Reference_Units)
  logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                " rows removed due to units listed in the Result Text varying from Reference Units"))
  fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_3_Mismatch_units", output_file_ending))
  Labs <- Labs %>%
    filter(Result_Text_Units == "" | Result_Text_Units == Reference_Units) %>% select(-Result_Text_Units)
  
  # Clean data (5)
  #   Get rid of # in Abnormal Flag because it tells you nothing (* on the otherhand means out of range);
  #   Replace LL with L and HH with H; Find "flags" in text and add to Abnormal_Flag column if missing
  Labs <- Labs %>% mutate(Abnormal_Flag = gsub("#", "", Abnormal_Flag),
                          Abnormal_Flag = gsub("(L{1,})", "L", Abnormal_Flag),
                          Abnormal_Flag = gsub("(H{1,})", "H", Abnormal_Flag),
                          Abnormal_Flag = ifelse(grepl("Flag: *H(igh)* ", Result_Text, ignore.case = TRUE) &
                                                   Abnormal_Flag == "", "H", Abnormal_Flag),
                          Abnormal_Flag = ifelse(grepl("Flag: *L(ow)* ", Result_Text, ignore.case = TRUE) &
                                                   Abnormal_Flag == "", "L", Abnormal_Flag))
  
  # Clean data (6)
  #   Split units and split date time
  #   Create AM/PM Flag
  Labs <- Labs %>%
    separate(Reference_Units, c("Unit_Num", "Unit_Den", sep = "/"), remove = FALSE) %>% select(-"/") %>%
    extract(Seq_Date_Time, c("Seq_Date", "Seq_Time"),
            regex = "(\\d{2}/\\d{2}/\\d{4}) (\\d{2}:\\d{2})", remove = FALSE) %>%
    separate(Seq_Time, c("Seq_Hour", "Seq_Min", sep = ":"), remove = FALSE) %>% select(-":") %>%
    mutate(Seq_Date = ifelse(is.na(Seq_Time), Seq_Date_Time, Seq_Date),
           Seq_Date = mdy(Seq_Date),
           Seq_Hour = as.numeric(Seq_Hour),
           Seq_Min = as.numeric(Seq_Min),
           FLAG_AM_PM = ifelse(is.na(Seq_Time), NA, ifelse(Seq_Hour < 12, "AM", "PM")))
  
  # Clean data (7)
  #   If strict, remove NA times, remove AM if PM in name, remove PM if AM in name
  if (strict){
    Lab_abn <- Labs %>% filter(is.na(Seq_Time))
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed because no Seq_Time specifed in Seq_Date_Time unit"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_4a_Strict_NA_Times", output_file_ending))
    Labs <- Labs %>% filter(!is.na(Seq_Time))
    
    Lab_abn <- Labs %>% filter(grepl("AM", Group_Id), FLAG_AM_PM != "AM")
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed because AM test not performed in AM"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_4b_Strict_AM_Times", output_file_ending))
    Labs <- Labs %>% filter((grepl("AM", Group_Id) & FLAG_AM_PM == "AM") | !grepl("AM", Group_Id))
    
    Lab_abn <- Labs %>% filter(grepl("PM", Group_Id), FLAG_AM_PM != "PM")
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed because PM test not performed in PM"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_4c_Strict_PM_Times", output_file_ending))
    Labs <- Labs %>% filter((grepl("PM", Group_Id) & FLAG_AM_PM == "PM") | !grepl("PM", Group_Id))
    
    Lab_abn <- Labs %>% filter(grepl("12.?am", Group_Id), Seq_Time != "00:00")
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed because 12am test not performed in 12am"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_4d_Strict_12am_Times", output_file_ending))
    Labs <- Labs %>% filter((grepl("12.{0,1}am", Group_Id) & Seq_Time == "00:00") | !grepl("12.{0,1}am", Group_Id))
  }
  
  Group_Id_list <- Labs %>% group_by(Group_Id) %>% summarise() %>% pull(Group_Id)
  unit_values <- c(dl = 1e-1, ml = 1e-3, ug = 1e-6, ng = 1e-9, pg = 1e-12)
  if (skip_ACTH) {Group_Id_list <- grep("^(?!ACTH).*$", Group_Id_list, value = TRUE, perl = TRUE)}
  # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
  cushings_threshold <- c(ug_dl = 1e2, ng_dl = 1e5, ug_ml = 1, ng_ml = 1e3)
  
  if(create_cortisol_group){
    Cortisol_Group_Id_list <- ifelse(skip_ACTH,
                                     grep("^Cortisol((?!ACTH).)*$", Group_Id_list, value = TRUE, perl = TRUE),
                                     grep("Cortisol", Group_Id_list, value = TRUE))
  }
  
  for (Id in Group_Id_list){
    header <- gsub("\\)", "", gsub("_+", "_", gsub("( |\\(|/|,)", "_", Id)))
    if (strict) { header <- str_c(header, "_Strict")}
    Subgroup <- Labs %>% filter(Group_Id == Id) %>% group_by(EMPI)
    logdebug(Id)
    logdebug(Subgroup %>% group_by(Reference_Units) %>% summarise(n = n()))
    
    # Note any possible duplicates, but don't necessarily remove
    nSubjects <- Subgroup %>% group_by(EMPI) %>% summarise(count = n()) %>% pull(count) %>% length()
    logdebug(str_c("Number of Subjects: ", nSubjects))
    select_EMPIs <- Subgroup %>% arrange(EMPI) %>% group_by(EMPI, Seq_Date_Time) %>%
      summarise(Count = n()) %>% filter(Count > 1) %>% pull(EMPI) %>% unique()
    if (length(select_EMPIs)){
      path_duplicates <- str_c(path_lab_abn, "Duplicate_date_time_EMPIs/")
      if(!dir.exists(path_duplicates)) {dir.create(path_duplicates)}
      logwarn(str_c(length(select_EMPIs), " out of ", nSubjects, " subjects have exact date and time duplicates"))
      Lab_abn <- Subgroup %>% filter(EMPI %in% select_EMPIs) %>% add_count(EMPI, Seq_Date_Time) %>% filter(n > 1)
      fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_5_Duplicate_date_time_", header, output_file_ending))
      for (empi in select_EMPIs){
        fwrite(Lab_abn %>% filter(EMPI == empi), str_c(path_duplicates, header, "_", empi, output_file_ending))
      }
      rm(empi, path_duplicates)
    }
    rm(nSubjects, select_EMPIs)
    
    # Find the main reference
    # - priority to reference unit included in name over the majority number if different
    main_unit <- ifelse(grepl("^.*\\(([[:alpha:]]{2,3}/[[:alpha:]]{2})\\).*$", Id),
                        tolower(gsub("^.*\\(([[:alpha:]]{2,3}/[[:alpha:]]{2})\\).*$", "\\1", Id)),
                        Subgroup %>% group_by(Reference_Units) %>% summarise(count = n()) %>%
                          filter(count == max(count)) %>% pull(Reference_Units))
    c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
    
    # Figure out if any reference range information is given as next few steps require for cleaning
    ref_ranges_summary <- Subgroup %>% group_by(Reference_Range) %>% summarise() %>% pull()
    if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
      logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Subgroup), " entries"))
      max_range = Subgroup %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
    } else {
      option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
      option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
      max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                       Subgroup %>% filter(Reference_Units == main_unit) %>%
                                         group_by(Reference_Range) %>% summarise() %>% pull())), na.rm = TRUE)
      rm(option1, option2)
    }
    rm(ref_ranges_summary)
    logdebug(str_c("max_range: ", max_range))
    
    # Clean data (8):
    # (A) If units are missing and above the general scope, remove
    # (B) If units are missing but in the general scope, add main_unit
    Lab_abn <- Subgroup %>% filter(Reference_Units == "", Result > max_range & !(grepl("H", Abnormal_Flag)))
    if (nrow(Lab_abn) > 0){
      logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Subgroup), " rows with GroupId ", Id,
                    " have been removed due to lacking reference units, falling out of the max range,",
                    " and having no marker for being outside the range"))
      fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_6_Missing_units_and_out_of_range_", header, output_file_ending))
      Subgroup <- Subgroup %>% filter(Reference_Units != "" | Result <= max_range | grepl("H", Abnormal_Flag))
    }
    Subgroup <- Subgroup %>% mutate(Unit_Num = ifelse(Reference_Units == "", main_unit_num, Unit_Num),
                                    Unit_Den = ifelse(Reference_Units == "", main_unit_den, Unit_Den),
                                    Reference_Units = ifelse(Reference_Units == "", main_unit, Reference_Units))
    
    # Clean data (9):
    # (A) Convert all units to the same unit
    # (B) Remove if changed units are above range and have not been flagged high
    # (C) Remove if main units are above range and have not been flagged high
    # (D) Remove if values are above cushing threshold (if Cortisol test)
    Subgroup <- Subgroup %>%
      mutate(Result_update = Result,
             Result_update = ifelse(Unit_Num != main_unit_num,
                                    Result_update * unit_values[Unit_Num] / unit_values[main_unit_num],
                                    Result_update),
             Result_update = ifelse(Unit_Den != main_unit_den,
                                    Result_update * unit_values[main_unit_den] / unit_values[Unit_Den],
                                    Result_update),
             Reference_Units_update = main_unit)
    Lab_abn <- Subgroup %>% filter(Reference_Units != main_unit, Result_update > max_range, Abnormal_Flag == "")
    if (nrow(Lab_abn) > 0){
      logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Subgroup), " rows with GroupId ", Id,
                    " removed due to possible incorrect reference units: original units not ", main_unit))
      fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_7A_Possible_incorrect_units_", header, output_file_ending))
      Subgroup <- Subgroup %>% filter(Reference_Units == main_unit | Result_update <= max_range | Abnormal_Flag != "")
    }
    Lab_abn <- Subgroup %>% filter(Reference_Units == main_unit, Result > max_range, Abnormal_Flag == "")
    if (nrow(Lab_abn) > 0){
      logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Subgroup), " rows with GroupId ", Id,
                    " removed due to possible incorrect reference units: ", main_unit, " should have been a different unit"))
      fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_7B_Possible_incorrect_units_", header, output_file_ending))
      Subgroup <- Subgroup %>% filter(Reference_Units != main_unit | Result <= max_range | Abnormal_Flag != "")
    }
    if (grepl("Cortisol", Id)){
      # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
      Lab_abn <- Subgroup %>% filter(Result <= cushings_threshold[main_unit], Result_update > cushings_threshold[main_unit])
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Subgroup), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: original result below cushing threshold but update is not"))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_7C_Possible_incorrect_units_", header, output_file_ending))
        Subgroup <- Subgroup %>% filter(Result > cushings_threshold[main_unit] | Result_update <= cushings_threshold[main_unit])
      }
      # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
      Lab_abn <- Subgroup %>% filter(Result_update > cushings_threshold[main_unit])
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Subgroup), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: result is above cushing threshold"))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_7D_Possible_incorrect_units_", header, output_file_ending))
        Subgroup <- Subgroup %>% filter(Result_update <= cushings_threshold[main_unit])
      }
    }
    Subgroup <- Subgroup %>% mutate(Result = Result_update, Reference_Units = Reference_Units_update) %>%
      select(-c(Result_update, Reference_Units_update))
    rm(max_range, main_unit, main_unit_num, main_unit_den, Lab_abn)
    
    if(write_files){
      fwrite(Subgroup, str_c(output_file_header, header, output_file_ending))
    }
    if (create_cortisol_group){
      if (Id %in% Cortisol_Group_Id_list){
        if (exists("Cortisol_group")){
          Cortisol_group <- rbind(Cortisol_group, Subgroup)
        } else {
          Cortisol_group <- Subgroup
        }
      }
    }
    
    DF_to_fill <- Create_ACTH_Cortisol_DHEA_Output_Columns(Subgroup, header, DF_to_fill)
    rm(header, Subgroup)
  }
  rm(Id)
  
  if (create_cortisol_group){
    header <- "All_Cortisol"
    if (strict) { header <- str_c(header, "_Strict")}
    logdebug(str_c("Number of Subjects: ", Cortisol_group %>% group_by(EMPI) %>% summarise(count = n()) %>% pull(count) %>% length()))
    # Cortisol should all be the same unit (usually ug/dl) but if that is not the case, change them to whichever unit is the most common
    if (Cortisol_group %>% group_by(Reference_Units) %>% summarise() %>% pull(Reference_Units) %>% length() > 1){
      logdebug(Cortisol_group %>% group_by(Reference_Units) %>% summarise(n = n()))
      # Find the main reference
      main_unit <- Cortisol_group %>% group_by(Reference_Units) %>% summarise(count = n()) %>%
        filter(count == max(count)) %>% pull(Reference_Units)
      c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
      # Figure out if any reference range information is given as next few steps require for cleaning
      ref_ranges_summary <- Cortisol_group %>% group_by(Reference_Range) %>% summarise() %>% pull()
      if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
        logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Cortisol_group), " entries"))
        max_range = Cortisol_group %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
      } else {
        option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
        option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
        max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                         Cortisol_group %>% filter(Reference_Units == main_unit) %>%
                                           group_by(Reference_Range) %>% summarise() %>% pull())), na.rm = TRUE)
        rm(option1, option2)
      }
      rm(ref_ranges_summary)
      logdebug(str_c("max_range: ", max_range))
      # Clean data (10):
      # (A) Convert all units to the same unit
      # (B) Remove if changed units are above range and have not been flagged high
      # (C) Remove if main units are above range and have not been flagged high
      # (D) Remove if values are above cushing threshold (if Cortisol test)
      Cortisol_group <- Cortisol_group %>%
        mutate(Result_update = Result,
               Result_update = ifelse(Unit_Num != main_unit_num,
                                      Result_update * unit_values[Unit_Num] / unit_values[main_unit_num],
                                      Result_update),
               Result_update = ifelse(Unit_Den != main_unit_den,
                                      Result_update * unit_values[main_unit_den] / unit_values[Unit_Den],
                                      Result_update),
               Reference_Units_update = main_unit)
      Lab_abn <- Cortisol_group %>% filter(Reference_Units != main_unit, Result_update > max_range, Abnormal_Flag == "")
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Cortisol_group), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: original units not ", main_unit))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_8A_Possible_incorrect_units_", header, output_file_ending))
        Cortisol_group <- Cortisol_group %>% filter(Reference_Units == main_unit | Result_update <= max_range | Abnormal_Flag != "")
      }
      Lab_abn <- Cortisol_group %>% filter(Reference_Units == main_unit, Result > max_range, Abnormal_Flag == "")
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Cortisol_group), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: ", main_unit, " should have been a different unit"))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_8B_Possible_incorrect_units_", header, output_file_ending))
        Cortisol_group <- Cortisol_group %>% filter(Reference_Units != main_unit | Result <= max_range | Abnormal_Flag != "")
      }
      # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
      Lab_abn <- Cortisol_group %>% filter(Result <= cushings_threshold[main_unit], Result_update > cushings_threshold[main_unit])
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Cortisol_group), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: original result below cushing threshold but update is not"))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_8C_Possible_incorrect_units_", header, output_file_ending))
        Cortisol_group <- Cortisol_group %>% filter(Result > cushings_threshold[main_unit] | Result_update <= cushings_threshold[main_unit])
      }
      # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
      Lab_abn <- Cortisol_group %>% filter(Result_update > cushings_threshold[main_unit])
      if (nrow(Lab_abn) > 0){
        logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Cortisol_group), " rows with GroupId ", Id,
                      " removed due to possible incorrect reference units: result is above cushing threshold"))
        fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_8D_Possible_incorrect_units_", header, output_file_ending))
        Cortisol_group <- Cortisol_group %>% filter(Result_update <= cushings_threshold[main_unit])
      }
      Cortisol_group <- Cortisol_group %>% mutate(Result = Result_update, Reference_Units = Reference_Units_update) %>%
        select(-c(Unit_Num, Unit_Den, Result_update, Reference_Units_update))
      rm(max_range, main_unit, main_unit_num, main_unit_den, Lab_abn)
    }
    
    if(write_files){
      fwrite(Cortisol_group, str_c(output_file_header, header, output_file_ending))
    }
    DF_to_fill <- Create_ACTH_Cortisol_DHEA_Output_Columns(Cortisol_group, header, DF_to_fill)
    
    rm(Cortisol_group, header)
  }
  rm(Labs)
  return(DF_to_fill)
}

process6_IGE_IGG_labs <- function(DF_to_fill = All_merged,
                                  input_file_header = config$rpdr_file_header_july,
                                  input_file_ending = config$rpdr_file_ending,
                                  path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                                  output_file_ending = config$general_file_ending){
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  
  # Clean data (1):
  # (A): Change "Less than x" to "< x" and "Greater than x" to "> x"; Get rid of an extra spacing
  # (B): Only include values that are digits, decimal starts, < x, > x
  Lab_abn <- Labs %>% filter(!grepl("^ *(<|>|LESS THAN|GREATER THAN)* *(\\d|\\.)", toupper(Result)))
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " entries out of ", nrow(Labs), " removed due to missingness or corrupt result entries"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_1_Missingness_or_corrupt_result", output_file_ending))
  }
  Labs <- Labs %>%
    mutate(Result = gsub("Less than", "<", Result, ignore.case = TRUE),
           Result = gsub("Greater than", ">", Result, ignore.case = TRUE),
           Result = gsub(" ", "", Result)) %>%
    filter(grepl("^(<|>)*(\\d|\\.)", Result))
  
  # Clean data (2)
  #   Change <0 to 0; Change all <x to x-min(x)/10; Change all >x to x + 0.11111
  Labs <- Labs %>%
    mutate(Result = gsub("<0((\\.|0)*)$", "0", Result),
           Result = gsub("(.*(\\d|\\.)+).*", "\\1", Result),
           LessThanX = as.numeric(ifelse(grepl("<", Result), gsub("<((\\d|\\.)+)", "\\1", Result), NA)),
           GreaterThanX = as.numeric(ifelse(grepl(">", Result), gsub(">((\\d|\\.)+)", "\\1", Result), NA)),
           Result = as.numeric(Result),
           Result = ifelse(is.na(LessThanX), Result, LessThanX - min(LessThanX, na.rm = TRUE) / 10),
           Result = ifelse(is.na(GreaterThanX), Result, GreaterThanX + 1/9)) %>%
    select(-c(LessThanX, GreaterThanX))
  
  # Clean data (3)
  #   Remove duplicates
  Lab_abn <- Labs[(duplicated(Labs)),]
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " completely duplicated row(s) out of ", nrow(Labs), " removed"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_2_Duplicate_rows", output_file_ending))
  }
  Labs <- Labs %>% unique()
  
  # Clean data (4)
  # (A): Change all reference units to upper case; IU/ML with KU/L (they are equal);
  #      Find units listed in result text
  # (B): Remove mismatches between units in result text and reference units
  Labs <- Labs %>%
    mutate(Reference_Units = toupper(Reference_Units),
           Reference_Units = ifelse(grepl("IU/ML", Reference_Units), "KU/L", Reference_Units),
           Result_Text = gsub("iu/ml", "KU/L", Result_Text, ignore.case = TRUE),
           Result_Text = gsub("mg/dl", "MG/DL", Result_Text, ignore.case = TRUE),
           Result_Text_Units = ifelse(grepl(".*(KU/L|MG/DL).*", Result_Text, ignore.case = TRUE),
                                      gsub(".*(KU/L|MG/DL).*", "\\1", Result_Text, ignore.case = TRUE),
                                      ""),
           Result_Text_Units = toupper(Result_Text_Units))
  Lab_abn <- Labs %>% filter(Result_Text_Units != "" & Result_Text_Units != Reference_Units)
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed due to units listed in the Result Text varying from Reference Units"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_3_Mismatch_units", output_file_ending))
  }
  Labs <- Labs %>%
    filter(Result_Text_Units == "" | Result_Text_Units == Reference_Units) %>% select(-Result_Text_Units)
  
  # Clean data (5)
  #   Remove units out of group 
  Lab_abn <- Labs %>% filter((Group_Id == "IGE" & Reference_Units != "KU/L") | 
                               (Group_Id == "IGG" & Reference_Units != "MG/DL"))
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows removed due to unit differing from group"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_4_Mismatch_units", output_file_ending))
  }
  Labs <- Labs %>% filter((Group_Id == "IGE" & Reference_Units == "KU/L") | 
                            (Group_Id == "IGG" & Reference_Units == "MG/DL"))
  
  # Clean data (6)
  #   Adjust date/time
  Labs <- Labs %>% mutate(Seq_Date_Time = mdy_hm(Seq_Date_Time)) %>% arrange(EMPI, Seq_Date_Time)
  
  # Clean data (7)
  #   Remove additional result from time point
  empis_with_mult_times <- Labs %>% filter(EMPI %in% (Labs %>% group_by(EMPI, Seq_Date_Time, Group_Id) %>%
                                                        summarise(Count = n()) %>% filter(Count > 1) %>% group_by(EMPI) %>% summarise() %>% pull())) %>%
    arrange(EMPI, Group_Id, Seq_Date_Time)
  relevant_info <- empis_with_mult_times %>% select(EMPI, Seq_Date_Time, Group_Id)
  Lab_abn <- empis_with_mult_times[duplicated(relevant_info)| duplicated(relevant_info, fromLast = TRUE), ]
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " out of ", nrow(Labs),
                  " rows have two similar values at the same timepoint"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_5_Duplicate_times", output_file_ending))
  }
  rm(empis_with_mult_times, relevant_info)
  rm(Lab_abn)
  
  Group_Id_list = c("IGE", "IGG")
  for (ID in Group_Id_list){
    Group <- Labs %>% filter(Group_Id == ID)
    Output_Columns <- Group %>% group_by(EMPI, Seq_Date_Time, Reference_Units) %>%
      summarise(mean_duplicates = mean(Result)) %>%
      group_by(EMPI, Reference_Units) %>%
      summarise(!!as.symbol(str_c(ID, "_nTotalValues")) := n(),
                !!as.symbol(str_c(ID, "_Mean")) := mean(mean_duplicates),
                !!as.symbol(str_c(ID, "_SD")) := sd(mean_duplicates),
                !!as.symbol(str_c(ID, "_All_Seq_Date_Times")) := paste(Seq_Date_Time, collapse = ";"),
                !!as.symbol(str_c(ID, "_All_Results")) := paste(mean_duplicates, collapse = ";")) %>%
      rename(!!as.symbol(str_c(ID, "_Reference_Units")) := Reference_Units)
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  }
  
  rm(Labs)
  return(DF_to_fill)
}
process6_Covid_labs <- function(DF_to_fill = All_merged,
                                input_file_header = config$rpdr_file_header,
                                input_file_ending = config$rpdr_file_ending,
                                path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                                output_file_ending = config$general_file_ending){
  loginfo("Processing Covid19 labs...")
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  
  Labs <- Labs %>% filter(grepl("SARS", Group_Id))
  
  Labs <- Labs %>% mutate(Result_Updated = ifelse(grepl("2 RNA", Group_Id),
                                                  # RNA Tests
                                                  ifelse(grepl("NEGATIVE|(NOT |UN)DETECTED", Result),
                                                         "NEGATIVE",
                                                         ifelse(grepl("POSITIVE|DETECTED", Result),
                                                                "POSITIVE",
                                                                NA)),
                                                  # IGG Tests
                                                  ifelse(grepl("POSITIVE", Result),
                                                         "POSITIVE",
                                                         ifelse(grepl("^Po", Result_Text),
                                                                "POSITIVE",
                                                                ifelse(grepl("^Ne", Result_Text),
                                                                       "NEGATIVE",
                                                                       ifelse(grepl("^<", Result),
                                                                              "NEGATIVE",
                                                                              NA))))))
  
  Lab_abn <- Labs[(duplicated(Labs)),]
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " completely duplicated row(s) out of ", nrow(Labs), " removed"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_Duplicate_rows", output_file_ending))
  }
  Labs <- Labs %>% unique()
  
  Labs <- Labs %>% filter(!is.na(Result_Updated)) %>% group_by(EMPI) %>% arrange(Seq_Date_Time)
  Output_Columns <- Labs %>% summarise("Tested" = "Yes",
                                       "Number_Of_Total_Tests" = n(),
                                       "First_Test_Date" = first(Seq_Date_Time),
                                       "First_Test_Result" = first(Result_Updated),
                                       "All_Test_Dates" = paste(Seq_Date_Time, collapse = ";"),
                                       "Result_Order" = paste(Result_Updated, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Tested = ifelse(is.na(Tested), "No", Tested),
                                      Number_Of_Total_Tests = ifelse(is.na(Number_Of_Total_Tests),
                                                                     0, Number_Of_Total_Tests))
  Subgroup <- Labs %>% filter(grepl("POSITIVE", Result_Updated))
  Output_Columns <- Subgroup %>% summarise("Positive_Result" = "Yes",
                                           "Number_Of_Total_Positive_Tests" = n(),
                                           "First_Positive_Test_Date" = first(Seq_Date_Time),
                                           "All_Positive_Test_Dates" = paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Positive_Result = ifelse(is.na(Positive_Result), "No", Positive_Result),
                                      Number_Of_Total_Positive_Tests = ifelse(is.na(Number_Of_Total_Positive_Tests),
                                                                              0, Number_Of_Total_Positive_Tests))
  Subgroup <- Labs %>% filter(grepl("NEGATIVE", Result_Updated))
  Output_Columns <- Subgroup %>% summarise("Negative_Result" = "Yes",
                                           "Number_Of_Total_Negative_Tests" = n(),
                                           "First_Negative_Test_Date" = first(Seq_Date_Time),
                                           "All_Negative_Test_Dates" = paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Negative_Result = ifelse(is.na(Negative_Result), "No", Negative_Result),
                                      Number_Of_Total_Negative_Tests = ifelse(is.na(Number_Of_Total_Negative_Tests),
                                                                              0, Number_Of_Total_Negative_Tests))
  
  rm(Labs, Subgroup)
  return(DF_to_fill)
}

process6_Cholesterol_labs <- function(DF_to_fill = All_merged,
                                input_file_header = config$rpdr_file_header,
                                input_file_ending = config$rpdr_file_ending,
                                path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                                output_file_ending = config$general_file_ending,
                                Optimal_Intermediate = 200,
                                Intermediate_High = 240){
  loginfo("Processing Cholesterol labs...")
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  
  Labs <- Labs %>% filter(grepl("Cholesterol", Group_Id))
  
  # Clean data (1): Remove non-numeric values
  Lab_abn <- Labs %>% filter(!grepl("^\\d", Result))
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " entries out of ", nrow(Labs), " removed due to non-numeric values"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_1_Non-numeric", output_file_ending))
    Labs <- Labs %>% filter(grepl("^\\d", Result))
  }
  Labs <- Labs %>% mutate(Result = as.numeric(Result))
  
  # Clean data (2): Remove duplicates
  Lab_abn <- Labs[(duplicated(Labs)),]
  if (nrow(Lab_abn) > 0){
    logwarn(str_c(nrow(Lab_abn), " completely duplicated row(s) out of ", nrow(Labs), " removed"))
    fwrite(Lab_abn, str_c(path_lab_abn, "Abnormality_2_Duplicate_rows", output_file_ending))
    Labs <- Labs %>% unique()
  }
  
  Labs <- Labs %>%
    extract(Seq_Date_Time, c("Seq_Date", "Seq_Time"),
            regex = "(\\d{2}/\\d{2}/\\d{4}) (\\d{2}:\\d{2})", remove = FALSE) %>%
    separate(Seq_Time, c("Seq_Hour", "Seq_Min", sep = ":"), remove = FALSE) %>% select(-":") %>%
    mutate(Seq_Date = ifelse(is.na(Seq_Time), Seq_Date_Time, Seq_Date),
           Seq_Date = mdy(Seq_Date),
           Seq_Hour = as.numeric(Seq_Hour),
           Seq_Min = as.numeric(Seq_Min)) %>%
    arrange(Seq_Date, Seq_Time)
  
  loginfo(str_c("Using ", Optimal_Intermediate, " and ", Intermediate_High,
                " as thresholds between Optimal, Intermediate, and High"))
  Output_Columns <- Labs %>% group_by(EMPI) %>%
    summarise(Cholesterol = "Yes",
              Cholesterol_number_recorded = n(),
              Cholesterol_average = mean(Result),
              Cholesterol_min = min(Result, na.rm = TRUE),
              Cholesterol_max = max(Result, na.rm = TRUE),
              Cholesterol_median = median(Result, na.rm = TRUE),
              Cholesterol_all_values = paste(Result, collapse = ";"),
              Cholesterol_all_dates = paste(Seq_Date_Time, collapse = ";"),
              Cholesterol_probable_category = ifelse(Cholesterol_median < Optimal_Intermediate,
                                                     "Optimal",
                                                     ifelse(Cholesterol_median < Intermediate_High,
                                                            "Intermediate", "High")))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Cholesterol = ifelse(is.na(Cholesterol), "No", Cholesterol),
           Cholesterol_number_recorded = ifelse(is.na(Cholesterol_number_recorded),
                                                0, Cholesterol_number_recorded))
  Optimal <- Labs %>% filter(Result < Optimal_Intermediate)
  Output_Columns <- Optimal %>% group_by(EMPI) %>%
    summarise(Cholesterol_optimal = "Yes",
              Cholesterol_optimal_number_recorded = n(),
              Cholesterol_optimal_values = paste(Result, collapse = ";"),
              Cholesterol_optimal_dates = paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Cholesterol_optimal = ifelse(is.na(Cholesterol_optimal), "No", Cholesterol_optimal),
           Cholesterol_optimal_number_recorded = ifelse(is.na(Cholesterol_optimal_number_recorded),
                                                    0, Cholesterol_optimal_number_recorded))
  Intermediate <- Labs %>% filter(Result >= Optimal_Intermediate & Result < Intermediate_High)
  Output_Columns <- Intermediate %>% group_by(EMPI) %>%
    summarise(Cholesterol_intermediate = "Yes",
              Cholesterol_intermediate_number_recorded = n(),
              Cholesterol_intermediate_values = paste(Result, collapse = ";"),
              Cholesterol_intermediate_dates = paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Cholesterol_intermediate = ifelse(is.na(Cholesterol_intermediate), "No", Cholesterol_intermediate),
           Cholesterol_intermediate_number_recorded = ifelse(is.na(Cholesterol_intermediate_number_recorded),
                                                        0, Cholesterol_intermediate_number_recorded))
  High <- Labs %>% filter(Result >= Intermediate_High)
  Output_Columns <- High %>% group_by(EMPI) %>%
    summarise(Cholesterol_high = "Yes",
              Cholesterol_high_number_recorded = n(),
              Cholesterol_high_values = paste(Result, collapse = ";"),
              Cholesterol_high_dates = paste(Seq_Date_Time, collapse = ";"))
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Cholesterol_high = ifelse(is.na(Cholesterol_high), "No", Cholesterol_high),
           Cholesterol_high_number_recorded = ifelse(is.na(Cholesterol_high_number_recorded),
                                                             0, Cholesterol_high_number_recorded))
  rm(Optimal, Intermediate, High, Labs, Output_Columns)
  return(DF_to_fill)
}
####################################################################################################
##################################### Health History functions #####################################
####################################################################################################
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

####################################################################################################
####################################### Encounters functions #######################################
####################################################################################################
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

####################################################################################################
########################################## Summary Stats  ##########################################
####################################################################################################
create_summary_stats_files <- function(Cleaned_DF = All_merged,
                                       summary_stats_dir = config$summary_stats_dir,
                                       additional_header = "",
                                       TS = timestamp,
                                       file_ending = config$general_file_ending,
                                       use.TS = TRUE,
                                       return.vital_status = TRUE,
                                       return.age = TRUE,
                                       return.race = TRUE,
                                       return.gender = TRUE,
                                       return.bmi = TRUE,
                                       return.cholesterol = TRUE,
                                       return.flu = TRUE){
  loginfo("Generate summaries...")
  stats_file_ending = ifelse(use.TS, str_c("_Stats_", TS, file_ending), str_c("_Stats", file_ending))
  stats_file_input = str_c(summary_stats_dir, additional_header)
  if (return.age){
    Cleaned_DF <- Cleaned_DF %>%
      mutate(Living = ifelse(grepl("Not reported as deceased", Vital_status), "Yes", "No"),
             Age_Range = ifelse(Age < 10,
                                "0-9",
                                ifelse(Age < 20,
                                       "10-19",
                                       ifelse(Age < 30,
                                              "20-29",
                                              ifelse(Age < 40,
                                                     "30-39",
                                                     ifelse(Age < 50,
                                                            "40-49",
                                                            ifelse(Age < 60,
                                                                   "50-59",
                                                                   ifelse(Age < 70,
                                                                          "60-69",
                                                                          ifelse(Age < 80,
                                                                                 "70-79",
                                                                                 ifelse(Age < 90,
                                                                                        "80-89",
                                                                                        "90+"))))))))))
  }
  if (!dir.exists(summary_stats_dir)) {dir.create(summary_stats_dir)}
  Alive <- Cleaned_DF %>% filter(grepl("Not reported as deceased", Vital_status))
  Deceased <- Cleaned_DF %>% filter(!grepl("Not reported as deceased", Vital_status))
  if (return.vital_status){
    loginfo("Creating vital status table...")
    fwrite(Cleaned_DF %>% group_by(Living) %>% summarise (Count_All = n(),
                                                          Percent_All = n()/nrow(Cleaned_DF)*100),
           str_c(stats_file_input, "Vital_Status", stats_file_ending))
  }
  if (return.gender){
    loginfo("Creating gender table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Gender) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Gender) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Gender) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Gender", stats_file_ending))
  }
  if (return.age){
    loginfo("Creating age table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Age_Range) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Age_Range) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Age_Range) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Age", stats_file_ending)) 
  }
  if (return.race){
    loginfo("Creating race table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Race) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Race) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Race) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Race", stats_file_ending))
  }
  if (return.bmi){
    loginfo("Creating bmi table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(BMI_probable_category) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(BMI_probable_category) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(BMI_probable_category) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "BMI", stats_file_ending))
  }
  if (return.cholesterol){
    loginfo("Creating cholesterol table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Cholesterol_probable_category) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Cholesterol_probable_category) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Cholesterol_probable_category) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Cholesterol", stats_file_ending))
  }
  if (return.flu){
    loginfo("Creating flu table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Flu_vaccined) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Flu_vaccined) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Flu_vaccined) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Flu_Vaccination", stats_file_ending))
  }
}

####################################################################################################
####################################### Pulmonary functions  #######################################
####################################################################################################
process8_pulmonary <- function(DF_to_fill = All_merged,
                               pul_file = config$pul_rdata_file_name,
                               path_dia_abn = str_c(config$data_dir, "Pulmonary_abnormalities/"),
                               write_files = config$create_intermediates,
                               output_file_header = config$intermediate_files_dir,
                               output_file_ending = config$general_file_ending){
  DF_to_fill <- DF_to_fill[,1:10]
  load(pul_file)
  pft.id <- data.frame(pft.id)
  pul.id_update <- data.frame(Orig = matrix(unlist(pul.id), nrow=length(pul.id), byrow=T)) %>%
    mutate(Count = str_count(Orig, "|"))
}
