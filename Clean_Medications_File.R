require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)
library(openxlsx)

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
  Med_Map_df <- read.xlsx(str_c(general_path, "/Medication_Mapping.xlsx"))
  logdebug("Available folders and number of medication groups:")
  logdebug(t(Med_Map_df %>% group_by(Medication_Biobank_Folder) %>%
               summarise(Count = n(),
                         .groups = 'drop')))
  # Select relevant files
  Relevant_Meds <- data.frame(Medication_Name = character(),
                              Medication_Biobank_Folder = character(),
                              Medication_GWAS_Group = character(),
                              Medication_Pegasus_Group = character(),
                              ICS_or_ICS_LABA = character(),
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
  Medications_condensed <- Medications %>% group_by(Medication) %>%
    summarise(.groups = 'drop')
  Med_all_subgroups <- data.frame(Medication = character(),
                                  Medication_Name = character(),
                                  Medication_Biobank_Folder = character(),
                                  Medication_GWAS_Group = character(),
                                  Medication_Pegasus_Group = character(),
                                  ICS_or_ICS_LABA = character())
  for (mri in Med_Row_ids){
    Med_subgroup <- Medications_condensed %>% filter(grepl(Relevant_Meds[mri,]$Search_Term,
                                                           Medication,
                                                           ignore.case = Relevant_Meds[mri,]$Ignore.case,
                                                           perl = Relevant_Meds[mri,]$Perl)) %>%
      mutate(Medication_Name  = Relevant_Meds[mri,]$Medication_Name,
             Medication_Biobank_Folder = Relevant_Meds[mri,]$Medication_Biobank_Folder,
             Medication_GWAS_Group = Relevant_Meds[mri,]$Medication_GWAS_Group,
             Medication_Pegasus_Group = Relevant_Meds[mri,]$Medication_Pegasus_Group,
             ICS_or_ICS_LABA = Relevant_Meds[mri,]$ICS_or_ICS_LABA)
    Med_all_subgroups <- rbind(Med_all_subgroups, Med_subgroup)
    rm(Med_subgroup)
  }
  Medications <- right_join(Medications, Med_all_subgroups, by = "Medication")
  rm(Med_all_subgroups, mri, Medications_condensed, Relevant_Meds, Med_Row_ids)
  # Med Cleanup
  Medications <- Medications %>% mutate(Medication_Date = mdy(Medication_Date)) %>% arrange(EMPI, Medication_Date)
  if (Daily_Dose_Info){
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
  }  
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
    Output_Columns = Medications %>% group_by(EMPI) %>%
      select(EMPI, Medication_Name) %>% unique() %>%
      summarise(!!(as.symbol(Merged_Group_Header)) := "Yes",
                .groups = 'drop')
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
    if (Group_Info){
      Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Medication_Name) %>% unique() %>%
        summarise(!!(as.symbol(Group_Header)) := "Yes",
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      DF_to_fill <- DF_to_fill %>%
        mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
      loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", Grouping_Name))
      rm(Output_Columns)
    }
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
        Output_Columns <- Subgroup %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(Subgroup_Header)) := "Yes",
                    .groups = 'drop')
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        DF_to_fill <- DF_to_fill %>%
          mutate(!!(as.symbol(Subgroup_Header)) := ifelse(is.na(!!(as.symbol(Subgroup_Header))), "No", "Yes"))
        Output_Columns <- Subgroup %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_first_date"))) := first(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_last_date"))) := last(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                    .groups = 'drop')
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        Output_Columns <- Subgroup %>% group_by(EMPI, Medication) %>%
          summarise(Medication_Occurances = n(),
                    .groups = 'drop') %>%
          group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                    !!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                    .groups = 'drop')
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        if (Daily_Dose_Info){
          Output_Columns <- Subgroup %>% group_by(EMPI) %>%
            summarise(!!(as.symbol(str_c(Subgroup_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                      .groups = 'drop')
          DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        }
        loginfo(str_c(nrow(Output_Columns), " subjects were prescribed ", tolower(Med_Name), "."))
        logdebug(str_c("nlines ", Subgroup_Header, " after completion: ", nrow(Subgroup)))
        if(write_files){
          fwrite(Subgroup, str_c(output_file_header, Subgroup_Header, output_file_ending))
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
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication_Name) %>%
        summarise(Medication_Occurances = n(),
                  .groups = 'drop') %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication) %>%
        summarise(Medication_Occurances = n(),
                  .groups = 'drop') %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      if (Daily_Dose_Info){
        Output_Columns <- Group %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                    .groups = 'drop')
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
                !!(as.symbol(str_c(Merged_Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                .groups = 'drop')
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication_Name) %>%
      summarise(Medication_Occurances = n(),
                .groups = 'drop') %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                .groups = 'drop')
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication) %>%
      summarise(Medication_Occurances = n(),
                .groups = 'drop') %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                .groups = 'drop')
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    if (Daily_Dose_Info){
      Output_Columns <- Medications_distinct %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    }
    if (nebs){
      Output_Columns <- Medications_distinct %>% filter(grepl("[Nn]eb", Medication)) %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) := "Yes",
                  .groups = 'drop')
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

process_medications_date_cutoff <- function(DF_to_fill_cutoff,
                                            input_file_header_cutoff,
                                            input_file_ending_cutoff,
                                            path_med_abn_cutoff,
                                            Medications_Of_Interest_cutoff,
                                            Group_Column_cutoff,
                                            Individual_Info_cutoff,
                                            Group_Info_cutoff,
                                            Daily_Dose_Info_cutoff,
                                            write_files_cutoff,
                                            output_file_header_cutoff,
                                            output_file_ending_cutoff,
                                            merged_group_name_cutoff,
                                            nebs_cutoff,
                                            date_variable_cutoff,
                                            restrict_to_before_cutoff){
  if (missing(Medications_Of_Interest_cutoff)){
    logerror("No list of Medications were specified. Process stopped.")
    return(DF_to_fill_cutoff)
  }
  # Usually this is put in the function definition but because it gets used in compare, needs to be defined here
  if (missing(DF_to_fill_cutoff)) {DF_to_fill_cutoff = All_merged}
  if (missing(input_file_header_cutoff)) {input_file_header_cutoff = config$rpdr_file_header}
  if (missing(input_file_ending_cutoff)) {input_file_ending_cutoff = config$rpdr_file_ending}
  if (missing(path_med_abn_cutoff)) {path_med_abn_cutoff = str_c(config$data_dir, "Medication_abnormalities/")}
  if (missing(Group_Column_cutoff)) {Group_Column_cutoff = config$medication_group}
  if (missing(Individual_Info_cutoff)) {Individual_Info_cutoff = TRUE}
  if (missing(Group_Info_cutoff)) {Group_Info_cutoff = TRUE}
  if (missing(Daily_Dose_Info_cutoff)) {Daily_Dose_Info_cutoff = FALSE}
  if (missing(write_files_cutoff)) {write_files_cutoff = config$create_intermediates}
  if (missing(output_file_header_cutoff)) {output_file_header_cutoff = config$intermediate_files_dir}
  if (missing(output_file_ending_cutoff)) {output_file_ending_cutoff = config$general_file_ending}
  if (missing(nebs_cutoff)) {nebs_cutoff = FALSE}
  if (missing(restrict_to_before_cutoff)) {restrict_to_before_cutoff = FALSE}
  
  loginfo("Processing medications file...")
  Medications <- fread(str_c(input_file_header_cutoff, "Med", input_file_ending_cutoff))
  if (!dir.exists(path_med_abn_cutoff)) {dir.create(path_med_abn_cutoff)}
  path_med_abn_cutoff <- str_c(path_med_abn_cutoff, ifelse(restrict_to_before_cutoff,
                                             "Med_PreDate_",
                                             "Med_PostDate_"))
  if(write_files_cutoff){
    output_file_header_cutoff <- ifelse(restrict_to_before_cutoff,
                                 str_c(output_file_header_cutoff, "Med_PreDate_"),
                                 str_c(output_file_header_cutoff, "Med_PostDate_"))
  }
  
  loginfo("Applying Med_Map...")
  Med_Map_df <- fread(str_c(general_path, "/Medication_Mapping.txt"))
  logdebug("Available folders and number of medication groups:")
  logdebug(t(Med_Map_df %>% group_by(Medication_Biobank_Folder) %>%
               summarise(Count = n(),
                         .groups = 'drop')))
  # Select relevant files
  Relevant_Meds <- data.frame(Medication_Name = character(),
                              Medication_Biobank_Folder = character(),
                              Medication_GWAS_Group = character(),
                              Medication_Pegasus_Group = character(),
                              Search_Term = character(),
                              Ignore.case = logical(),
                              Perl = logical())
  for (Grouping_Name in names(Medications_Of_Interest_cutoff)){
    if (length(Medications_Of_Interest_cutoff[[Grouping_Name]]) == 0){
      Mini_Med_Map <- Med_Map_df %>% filter(!!(as.symbol(Group_Column_cutoff)) %in% Grouping_Name)
      Medications_Of_Interest_cutoff[[Grouping_Name]] = Mini_Med_Map %>% pull(Medication_Name)
    } else {
      Mini_Med_Map <- Med_Map_df %>% filter(!!(as.symbol(Group_Column_cutoff)) %in% Grouping_Name,
                                            Medication_Name %in% Medications_Of_Interest_cutoff[[Grouping_Name]])
    }
    if (nrow(Mini_Med_Map) == 0){
      logerror("Med grouping does not exist. Please try again.")
      return(DF_to_fill_cutoff)
    }
    Relevant_Meds <- rbind(Relevant_Meds, Mini_Med_Map)
  }
  rm(Mini_Med_Map, Med_Map_df)
  Med_Row_ids <- 1:nrow(Relevant_Meds)
  Medications_condensed <- Medications %>% group_by(Medication) %>%
    summarise(.groups = 'drop')
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
  loginfo("Restricting prescription by cutoff...")
  # Restrict ids to only the ones in the reduced id list
  Medications <- Medications %>% filter(EMPI %in% DF_to_fill_cutoff$EMPI)
  Medications <- Medications %>% mutate(Medication_Date = mdy(Medication_Date)) %>%
    arrange(EMPI, Medication_Date)
  if (!sum(str_count(DF_to_fill_cutoff %>% pull(date_variable_cutoff)) == 10) == nrow(DF_to_fill_cutoff)){
    EMPI_Date_Limit <- DF_to_fill_cutoff %>% select(EMPI, date_variable_cutoff) %>%
      rename(Cutoff_Date_Time = date_variable_cutoff) %>%
      extract(Cutoff_Date_Time, c("Cutoff_Date", "Cutoff_Time"),
              regex = "(\\d{4}-\\d{2}-\\d{2}) (\\d{2}:\\d{2})", remove = FALSE) %>%
      mutate(Cutoff_Date = ifelse(is.na(Cutoff_Time), Cutoff_Date_Time, Cutoff_Date)) %>% select(EMPI, Cutoff_Date)
  } else {
    EMPI_Date_Limit <- DF_to_fill_cutoff %>% select(EMPI, date_variable_cutoff) %>%
      rename(Cutoff_Date = date_variable_cutoff)
  }
  # Restrict Medications to only before or only after
  Medications <- left_join(Medications, EMPI_Date_Limit)
  time_str = ifelse(restrict_to_before_cutoff, "before", "after")
  loginfo(str_c("Restricting medications to ", time_str, " cutoff dates... "))

  Original_row_length <- nrow(Medications)
  if(restrict_to_before_cutoff){
    Medications <- Medications %>% filter(Medication_Date <= Cutoff_Date)
  } else {
    Medications <- Medications %>% filter(Medication_Date >= Cutoff_Date)
  }
  New_row_length <- nrow(Medications)
  logdebug(str_c(Original_row_length, " medications reduced to ", New_row_length, " medications."))
  rm(time_str, Original_row_length, New_row_length)
  
  if (Daily_Dose_Info_cutoff){
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
  }  
  if ("Medication_Date_Detail" %in% names(Medications)){
    Med_abn <- Medications %>% filter(Medication_Date_Detail == "Removed")
    if (nrow(Med_abn) > 0){
      logwarn(str_c(nrow(Med_abn), " row(s) out of ", nrow(Medications), " have been removed due having the flag 'Removed'."))
      fwrite(Med_abn, str_c(path_med_abn_cutoff, "Abnormality_1_Removed_Flag_All_Medications", output_file_ending_cutoff))
      Medications <- Medications %>% filter(Medication_Date_Detail != "Removed")
    }
  }
  loginfo("Creating medication output columns...")
  # Get the "Any exist" first to lower the search group/increase speed later
  if (!missing(merged_group_name_cutoff)){
    Merged_Group_Header = str_c(ifelse(restrict_to_before_cutoff, "Med_PreDate_", "Med_PostDate_"),
                                "Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", merged_group_name_cutoff)))
    Output_Columns = Medications %>% group_by(EMPI) %>%
      select(EMPI, Medication_Name) %>% unique() %>%
      summarise(!!(as.symbol(Merged_Group_Header)) := "Yes",
                .groups = 'drop')
    DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    DF_to_fill_cutoff <- DF_to_fill_cutoff %>%
      mutate(!!(as.symbol(Merged_Group_Header)) := ifelse(is.na(!!(as.symbol(Merged_Group_Header))),
                                                          "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", merged_group_name_cutoff,
                  ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction"))
  }
  for (Grouping_Name in names(Medications_Of_Interest_cutoff)){
    Group_Header = str_c(ifelse(restrict_to_before_cutoff, "Med_PreDate_", "Med_PostDate_"),
                         "Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)))
    Group <- Medications %>%
      filter(grepl(Grouping_Name, !!(as.symbol(Group_Column_cutoff))),
             Medication_Name %in% Medications_Of_Interest_cutoff[[Grouping_Name]])
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Medication_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes",
                .groups = 'drop')
    if (Group_Info_cutoff){
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      DF_to_fill_cutoff <- DF_to_fill_cutoff %>%
        mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
    }
    loginfo(str_c(nrow(Output_Columns), " subjects were prescribed any ", Grouping_Name,
                  ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction"))
    rm(Output_Columns)
    if(Individual_Info_cutoff){
      # Look for the individual prescriptions
      for (Med_Name in Medications_Of_Interest_cutoff[[Grouping_Name]]){
        Subgroup_Header <- str_c(ifelse(restrict_to_before_cutoff, "Med_PreDate_", "Med_PostDate_"),
                                 gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)),
                                 "_", gsub("/| ", "_", Med_Name))
        Subgroup <- Group %>% filter(Medication_Name == Med_Name) %>% group_by(EMPI)
        Med_abn <- Subgroup[duplicated(Subgroup) | duplicated(Subgroup, fromLast = TRUE),]
        if (nrow(Med_abn) > 0){
          logwarn(str_c(nrow(Med_abn), " completely duplicated row(s) out of ", nrow(Subgroup), " found. Duplicates removed."))
          fwrite(Med_abn, str_c(path_med_abn_cutoff, "Abnormality_2_Duplicate_rows_", Subgroup_Header, output_file_ending_cutoff))
          Subgroup <- Subgroup %>% unique()
        }
        Subgroup2 <- Subgroup %>% select(EMPI, Medication_Date, Medication, Hospital)
        Med_abn <- Subgroup[duplicated(Subgroup2) | duplicated(Subgroup2, fromLast = TRUE),]
        if (nrow(Med_abn) > 0){
          logwarn(str_c(nrow(Med_abn), " partially duplicated row(s) out of ", nrow(Subgroup), " found. ",
                        sum(duplicated(Subgroup2, fromLast = TRUE)), " rows removed."))
          fwrite(Med_abn, str_c(path_med_abn_cutoff, "Abnormality_3_Duplicate_rows_", Subgroup_Header, output_file_ending_cutoff))
          Subgroup <- distinct(Subgroup, EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup2)
        # Start writing output
        Output_Columns <- Subgroup %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(Subgroup_Header)) := "Yes",
                    .groups = 'drop')
        DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
        DF_to_fill_cutoff <- DF_to_fill_cutoff %>%
          mutate(!!(as.symbol(Subgroup_Header)) := ifelse(is.na(!!(as.symbol(Subgroup_Header))), "No", "Yes"))
        Output_Columns <- Subgroup %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_first_date"))) := first(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_last_date"))) := last(Medication_Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                    .groups = 'drop')
        DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
        Output_Columns <- Subgroup %>% group_by(EMPI, Medication) %>%
          summarise(Medication_Occurances = n(),
                    .groups = 'drop') %>%
          group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                    !!(as.symbol(str_c(Subgroup_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                    .groups = 'drop')
        DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
        if (Daily_Dose_Info_cutoff){
          Output_Columns <- Subgroup %>% group_by(EMPI) %>%
            summarise(!!(as.symbol(str_c(Subgroup_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                      !!(as.symbol(str_c(Subgroup_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                      .groups = 'drop')
          DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
        }
        loginfo(str_c(nrow(Output_Columns), " subjects were prescribed ", tolower(Med_Name),
                      ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction", "."))
        logdebug(str_c("nlines ", Subgroup_Header, " after completion: ", nrow(Subgroup)))
        if(write_files_cutoff){
          fwrite(Subgroup, str_c(output_file_header_cutoff, Subgroup_Header, output_file_ending_cutoff))
        }
        rm(Subgroup_Header)
        rm(Subgroup, Output_Columns, Med_abn)
      }
      rm(Med_Name)
    }
    if(Group_Info_cutoff){
      # Now that all the errors have been noted.. add more count/date columns
      Group <- Group %>% distinct(EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
      Output_Columns <- Group %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_total_dates"))) := n(),
                  !!(as.symbol(str_c(Group_Header, "_first_date"))) := first(Medication_Date),
                  !!(as.symbol(str_c(Group_Header, "_last_date"))) := last(Medication_Date),
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication_Name) %>%
        summarise(Medication_Occurances = n(),
                  .groups = 'drop') %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      Output_Columns <- Group %>% group_by(EMPI, Medication) %>%
        summarise(Medication_Occurances = n(),
                  .groups = 'drop') %>%
        group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                  !!(as.symbol(str_c(Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      if (Daily_Dose_Info_cutoff){
        Output_Columns <- Group %>% group_by(EMPI) %>%
          summarise(!!(as.symbol(str_c(Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                    !!(as.symbol(str_c(Group_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                    .groups = 'drop')
        DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      }
      # Rearrange so "Any" columns are together for easier viewing/understanding
      DF_to_fill_cutoff <- DF_to_fill_cutoff %>% select(EMPI:Group_Header, starts_with(Group_Header), everything())
      rm(Output_Columns)
      if(write_files_cutoff){
        fwrite(Group, str_c(output_file_header_cutoff, Group_Header, output_file_ending_cutoff))
      }
    }
    rm(Group_Header)
    rm(Group)
  }
  if (!missing(merged_group_name_cutoff)){
    Medications_distinct <- Medications %>% distinct(EMPI, Medication_Date, Medication, Hospital, .keep_all = TRUE)
    Merged_Group_Header = str_c(ifelse(restrict_to_before_cutoff, "Med_PreDate_", "Med_PostDate_"),
                                "Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", merged_group_name_cutoff)))
    Output_Columns <- Medications_distinct %>% arrange(Medication_Date) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_total_dates"))) := n(),
                !!(as.symbol(str_c(Merged_Group_Header, "_first_date"))) := first(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_last_date"))) := last(Medication_Date),
                !!(as.symbol(str_c(Merged_Group_Header, "_dates"))) := paste(Medication_Date, collapse = ";"),
                .groups = 'drop')
    DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication_Name) %>%
      summarise(Medication_Occurances = n(),
                .groups = 'drop') %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_type"))) := paste(Medication_Name[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                .groups = 'drop')
    DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    Output_Columns <- Medications_distinct %>% group_by(EMPI, Medication) %>%
      summarise(Medication_Occurances = n(),
                .groups = 'drop') %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription_total"))) := max(Medication_Occurances),
                !!(as.symbol(str_c(Merged_Group_Header, "_most_common_prescription"))) := paste(Medication[which(Medication_Occurances == max(Medication_Occurances))], collapse = ";"),
                .groups = 'drop')
    DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    if (Daily_Dose_Info_cutoff){
      Output_Columns <- Medications_distinct %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_all_medications"))) := paste(Medication, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose_MCG"))) := paste(DAILY_DOSE_MCG, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_puffs"))) := paste(DAILY_DOSE_PUFF, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_daily_dose"))) := paste(DAILY_DOSE, collapse = "|"),
                  !!(as.symbol(str_c(Merged_Group_Header, "_notes"))) := paste(NOTES, collapse= "|"),
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    }
    if (nebs_cutoff){
      Output_Columns <- Medications_distinct %>% filter(grepl("[Nn]eb", Medication)) %>% group_by(EMPI) %>%
        summarise(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) := "Yes",
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      DF_to_fill_cutoff <- DF_to_fill_cutoff %>%
        mutate(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer"))) :=
                 ifelse(is.na(!!(as.symbol(str_c(Merged_Group_Header, "_Nebulizer")))),
                        "No", "Yes"))
    }
    # Rearrange so "Any" columns are together for easier viewing/understanding
    DF_to_fill_cutoff <- DF_to_fill_cutoff %>% select(EMPI:Merged_Group_Header, starts_with(Merged_Group_Header), everything())
    rm(Output_Columns)
    if(write_files_cutoff){
      fwrite(Medications_distinct, str_c(output_file_header_cutoff, Merged_Group_Header, output_file_ending_cutoff))
    }
    rm(Medications_distinct)
  }
  return(DF_to_fill_cutoff)
}

process_medications_date_compare_cutoff <- function(DF_to_fill = ICS_Tested,
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
                                                    nebs = FALSE,
                                                    cutoff_variable){
  DF_to_fill <- process_medications_date_cutoff(DF_to_fill_cutoff = DF_to_fill,
                                                input_file_header_cutoff = input_file_header,
                                                input_file_ending_cutoff = input_file_ending,
                                                path_med_abn_cutoff = path_med_abn,
                                                Medications_Of_Interest_cutoff = Medications_Of_Interest,
                                                Group_Column_cutoff = Group_Column,
                                                Individual_Info_cutoff = Individual_Info,
                                                Group_Info_cutoff = Group_Info,
                                                Daily_Dose_Info_cutoff = Daily_Dose_Info,
                                                write_files_cutoff = write_files,
                                                output_file_header_cutoff = output_file_header,
                                                output_file_ending_cutoff = output_file_ending,
                                                merged_group_name_cutoff = merged_group_name,
                                                nebs_cutoff = nebs,
                                                date_variable_cutoff = cutoff_variable,
                                                restrict_to_before_cutoff = TRUE)
  DF_to_fill <- process_medications_date_cutoff(DF_to_fill_cutoff = DF_to_fill,
                                                input_file_header_cutoff = input_file_header,
                                                input_file_ending_cutoff = input_file_ending,
                                                path_med_abn_cutoff = path_med_abn,
                                                Medications_Of_Interest_cutoff = Medications_Of_Interest,
                                                Group_Column_cutoff = Group_Column,
                                                Individual_Info_cutoff = Individual_Info,
                                                Group_Info_cutoff = Group_Info,
                                                Daily_Dose_Info_cutoff = Daily_Dose_Info,
                                                write_files_cutoff = write_files,
                                                output_file_header_cutoff = output_file_header,
                                                output_file_ending_cutoff = output_file_ending,
                                                merged_group_name_cutoff = merged_group_name,
                                                nebs_cutoff = nebs,
                                                date_variable_cutoff = cutoff_variable,
                                                restrict_to_before_cutoff = FALSE)
  Pre_string <- "Med_PreDate_"
  Post_string <- "Med_PostDate_"
  if (!missing(merged_group_name)){
    Merged_Group_Header <- str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", merged_group_name)))
    Compare_Header <- str_c("Compare_", Merged_Group_Header)
    Pre_Header <- str_c(Pre_string, Merged_Group_Header)
    Post_Header <- str_c(Post_string, Merged_Group_Header)
    DF_to_fill <- DF_to_fill %>%
      mutate(!!as.symbol(Compare_Header) := 
               ifelse(!!as.symbol(Pre_Header) == "Yes",
                      ifelse(!!as.symbol(Post_Header) == "Yes",
                             "Yes - Both ranges",
                             "Yes - Predates cutoff"),
                      ifelse(!!as.symbol(Post_Header) == "Yes",
                             "Yes - Postdates cutoff",
                             "No")),
             !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) := 
               ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type"))),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type"))),
                             NA,
                             "Post only"),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type"))),
                             "Pre only",
                             ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type")) ==
                                      !!as.symbol(str_c(Post_Header, "_most_common_prescription_type")),
                                    "Same",
                                    "Different"))),
             !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type_total")) :=
               ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total"))),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total"))),
                             "No Change",
                             "Increased"),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total"))),
                             "Decreased",
                             ifelse(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total")) >
                                      !!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total")),
                                    "Increased",
                                    ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total")) >
                                             !!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total")),
                                           "Decreased",
                                           "No Change")))),
             !!as.symbol(str_c(Compare_Header, "_most_common_prescription")) := 
               ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription"))),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                             NA,
                             "Post only"),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                             "Pre only",
                             ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription")) ==
                                      !!as.symbol(str_c(Post_Header, "_most_common_prescription")),
                                    "Same",
                                    "Different"))),
             !!as.symbol(str_c(Compare_Header, "_most_common_prescription_total")) :=
               ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total"))),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                             "No Change",
                             "Increased"),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                             "Decreased",
                             ifelse(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total")) >
                                      !!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")),
                                    "Increased",
                                    ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")) >
                                             !!as.symbol(str_c(Post_Header, "_most_common_prescription_total")),
                                           "Decreased",
                                           "No Change")))))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                  " subjects were only prescribed medications before their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                  " subjects were only prescribed medications after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                  " subjects were prescribed medications before and after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) == "Different")),
                  " subjects had a most common prescription type change before and after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) == "Same")),
                  " subjects most common prescription type stayed the same before and after their test"))
  }
  for (Grouping_Name in names(Medications_Of_Interest)){
    Group_Header = str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)))
    if (Group_Info){
      Compare_Header <- str_c("Compare_", Group_Header)
      Pre_Header <- str_c(Pre_string, Group_Header)
      Post_Header <- str_c(Post_string, Group_Header)
      DF_to_fill <- DF_to_fill %>%
        mutate(!!as.symbol(Compare_Header) := 
                 ifelse(!!as.symbol(Pre_Header) == "Yes",
                        ifelse(!!as.symbol(Post_Header) == "Yes",
                               "Yes - Both ranges",
                               "Yes - Predates cutoff"),
                        ifelse(!!as.symbol(Post_Header) == "Yes",
                               "Yes - Postdates cutoff",
                               "No")),
               !!as.symbol(str_c(Compare_Header, "_total_dates")) :=
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_total_dates"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_dates"))),
                               "No Change",
                               "Increased"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_dates"))),
                               "Decreased",
                               ifelse(!!as.symbol(str_c(Post_Header, "_total_dates")) >
                                        !!as.symbol(str_c(Pre_Header, "_total_dates")),
                                      "Increased",
                                      ifelse(!!as.symbol(str_c(Pre_Header, "_total_dates")) >
                                               !!as.symbol(str_c(Post_Header, "_total_dates")),
                                             "Decreased",
                                             "No Change")))),
               !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) := 
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type"))),
                               NA,
                               "Post only"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type"))),
                               "Pre only",
                               ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type")) ==
                                        !!as.symbol(str_c(Post_Header, "_most_common_prescription_type")),
                                      "Same",
                                      "Different"))),
               !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type_total")) :=
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total"))),
                               "No Change",
                               "Increased"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total"))),
                               "Decreased",
                               ifelse(!!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total")) >
                                        !!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total")),
                                      "Increased",
                                      ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_type_total")) >
                                               !!as.symbol(str_c(Post_Header, "_most_common_prescription_type_total")),
                                             "Decreased",
                                             "No Change")))),
               !!as.symbol(str_c(Compare_Header, "_most_common_prescription")) := 
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                               NA,
                               "Post only"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                               "Pre only",
                               ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription")) ==
                                        !!as.symbol(str_c(Post_Header, "_most_common_prescription")),
                                      "Same",
                                      "Different"))),
               !!as.symbol(str_c(Compare_Header, "_most_common_prescription_total")) :=
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                               "No Change",
                               "Increased"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                               "Decreased",
                               ifelse(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total")) >
                                        !!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")),
                                      "Increased",
                                      ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")) >
                                               !!as.symbol(str_c(Post_Header, "_most_common_prescription_total")),
                                             "Decreased",
                                             "No Change")))))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                    " subjects were only prescribed ", Med_Name, " before their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                    " subjects were only prescribed ", Med_Name, " after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                    " subjects were prescribed ", Med_Name, " before and after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "Increased")),
                    " subjects had an increase in ", Med_Name, " prescriptions"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "Decreased")),
                    " subjects had an decrease in ", Med_Name, " prescriptions"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "No Change")),
                    " subjects had no change in the number of ", Med_Name, " prescriptions"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) == "Different")),
                    " subjects had a most common prescription type change before and after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_most_common_prescription_type")) == "Same")),
                    " subjects most common prescription type stayed the same before and after their test"))
    }
    if(Individual_Info){
      # Look for the individual prescriptions
      for (Med_Name in Medications_Of_Interest[[Grouping_Name]]){
        Subgroup_Header <- str_c(gsub("_{1,}", "_", gsub(" |,|-", "_", Grouping_Name)),
                                 "_", gsub("/| ", "_", Med_Name))
        Compare_Header <- str_c("Compare_", Subgroup_Header)
        Pre_Header <- str_c(Pre_string, Subgroup_Header)
        Post_Header <- str_c(Post_string, Subgroup_Header)
        DF_to_fill <- DF_to_fill %>%
          mutate(!!as.symbol(Compare_Header) := 
                   ifelse(!!as.symbol(Pre_Header) == "Yes",
                          ifelse(!!as.symbol(Post_Header) == "Yes",
                                 "Yes - Both ranges",
                                 "Yes - Predates cutoff"),
                          ifelse(!!as.symbol(Post_Header) == "Yes",
                                 "Yes - Postdates cutoff",
                                 "No")),
                 !!as.symbol(str_c(Compare_Header, "_total_dates")) :=
                   ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_total_dates"))),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_dates"))),
                                 "No Change",
                                 "Increased"),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_dates"))),
                                 "Decreased",
                                 ifelse(!!as.symbol(str_c(Post_Header, "_total_dates")) > !!as.symbol(str_c(Pre_Header, "_total_dates")),
                                        "Increased",
                                        ifelse(!!as.symbol(str_c(Pre_Header, "_total_dates")) > !!as.symbol(str_c(Post_Header, "_total_dates")),
                                               "Decreased",
                                               "No Change")))),
                 !!as.symbol(str_c(Compare_Header, "_most_common_prescription")) := 
                   ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription"))),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                                 NA,
                                 "Post only"),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription"))),
                                 "Pre only",
                                 ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription")) ==
                                          !!as.symbol(str_c(Post_Header, "_most_common_prescription")),
                                        "Same",
                                        "Different"))),
                 !!as.symbol(str_c(Compare_Header, "_most_common_prescription_total")) :=
                   ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total"))),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                                 "No Change",
                                 "Increased"),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total"))),
                                 "Decreased",
                                 ifelse(!!as.symbol(str_c(Post_Header, "_most_common_prescription_total")) > !!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")),
                                        "Increased",
                                        ifelse(!!as.symbol(str_c(Pre_Header, "_most_common_prescription_total")) > !!as.symbol(str_c(Post_Header, "_most_common_prescription_total")),
                                               "Decreased",
                                               "No Change")))))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                      " subjects were only prescribed ", Med_Name, " before their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                      " subjects were only prescribed ", Med_Name, " after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                      " subjects were prescribed ", Med_Name, " before and after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "Increased")),
                      " subjects had an increase in ", Med_Name, " prescriptions"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "Decreased")),
                      " subjects had an decrease in ", Med_Name, " prescriptions"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_dates")) == "No Change")),
                      " subjects had no change in the number of ", Med_Name, " prescriptions"))
        rm(Subgroup_Header)
      }
      rm(Med_Name)
    }
    rm(Group_Header)
  }
  return(DF_to_fill)
}