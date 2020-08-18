require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_medications <- function(DF_to_fill = All_merged,
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
  if(write_files){
    output_file_header <- str_c(output_file_header, "Med_")
  }
  
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