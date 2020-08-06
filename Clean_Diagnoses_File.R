require(data.table) # fread, fwrite
require(dplyr) # filter, mutate, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

process_diagnoses <- function(DF_to_fill = All_merged,
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