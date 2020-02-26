require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(data.table) # fread, fwrite
require(tidyverse)
require(logging)
require(readr)
require(zeallot) # %<-% (multiple variable assignment)
require(sqldf)

####################################################################################################
######################################## General functions  ########################################
####################################################################################################
process1_biobankIDs <- function(input_file_header = config$rpdr_file_header,
                               input_file_ending = config$rpdr_file_ending){
  loginfo("Processing biobank ids file... ")
  BiobankIDs <- data.table(fread(str_c(input_file_header, "Bib", input_file_ending)))
  BiobankIDs <- BiobankIDs %>% select(Subject_Id, EMPI) %>% rename(Biobank_Subject_ID = Subject_Id)
  loginfo(str_c(nrow(BiobankIDs), " subjects processed"))
  return(BiobankIDs)
}

process2_demographics <- function(DF_to_fill = All_merged,
                                 input_file_header = config$rpdr_file_header,
                                 input_file_ending = config$rpdr_file_ending){
  loginfo("Processing demographics file...")
  Demographics <- data.table(fread(str_c(input_file_header, "Dem", input_file_ending)))
  Demographics <- Demographics %>% select(EMPI, Gender, Date_of_Birth, Age, Date_Of_Death, Race)
  DF_to_fill <- merge(DF_to_fill, Demographics, by = "EMPI")
  loginfo("Gender, date of birth, age, date of death, and race information have been added")
  rm(Demographics)
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
    mutate_at(vars(contains("Asthma")), ~ifelse(.x == "Yes", 1, 0))
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
    select(Biobank_Subject_ID, Asthma, NonAsthma, COPD, NonCOPD)
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and COPD identifiers have been added")
  rm(Deidentified)
  return(DF_to_fill)
}


process_BMI <- function(BMI_df = Deidentified){
  loginfo("Processing BMI...")
  BMI_df <- BMI_df %>% mutate_at(vars(contains("List")), ~gsub("\\[((\\d|\\.)+)\\]", "\\1", .x))
  BMI_df <- BMI_df %>% mutate_at(vars(contains("List")), ~gsub(" ", ";", .x))
  BMI_df <- BMI_df %>% mutate_at(vars(contains("List")), ~gsub("^$", NA, .x))
  return(BMI_df)
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
    select(-(starts_with("BMI")))
  Deidentified <- process_BMI(Deidentified)
  Deidentified <- process_dates(Deidentified)
  DF_to_fill <- left_join(DF_to_fill, Deidentified, by = "Biobank_Subject_ID")
  loginfo("Asthma and NonAsthma, BMI, and Spirometry identifiers have been added")
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
    if(write_files){
      fwrite(Group, str_c(output_file_header, Group_Header, output_file_ending))
    }
    rm(Group_Header)
    rm(Group)
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
                                 merged_group_name){
  if (missing(Medications_Of_Interest)){
    logerror("No list of Medications were specified. Process stopped.")
    return(DF_to_fill)
  }
  loginfo("Processing medications file...")
  Medications <- fread(str_c(input_file_header, "Med", input_file_ending))
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
    Output_Columns = Medications_distinct %>% group_by(EMPI) %>% select(EMPI, Medication_Name) %>%
      unique() %>% summarise(!!(as.symbol(Merged_Group_Header)) := "Yes")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(!!(as.symbol(Merged_Group_Header)) := ifelse(is.na(!!(as.symbol(Merged_Group_Header))),
                                                          "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", merged_group_name))
    Output_Columns <- Medications_distinct %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_total_dates"))) := n(),
                !!(as.symbol(str_c(Merged_Group_Header, "_first_date"))) := first(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_last_date"))) := last(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"))
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% filter(grepl("[Nn]eb", Medication)) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) := "Yes")
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>%
      mutate(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) := ifelse(is.na(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer")))),
                                                          "No", "Yes"))
      
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