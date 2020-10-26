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
                              Exact = FALSE,
                              Individual_Info = TRUE,
                              Group_Info = TRUE,
                              Merge_Group_Info_Name,
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
  # store relevant
  if(!missing(Merge_Group_Info_Name)){
    Relevant_rows <- NULL
    Original_Columns <- names(DF_to_fill)
  }
  if(write_files){
    output_file_header <- str_c(output_file_header, "Dia_")
  }
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  # Get the "Any exist" first to lower the search group/increase speed later
  for (Grouping_Name in names(Diagnoses_Of_Interest)){
    Group_Header = str_c("Any_", gsub(" ", "_", Grouping_Name))
    Group <- Diagnoses %>%
      filter(grepl(str_c(Diagnoses_Of_Interest[[Grouping_Name]], collapse = "|"), Diagnosis_Name))
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Diagnosis_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes",
                .groups = 'drop')
    if (Group_Info){
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      DF_to_fill <- DF_to_fill %>%
        mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
    }
    loginfo(str_c(nrow(Output_Columns), " subjects have any ", Grouping_Name, " diagnosis"))
    rm(Output_Columns)
    
    if (Individual_Info){
      if(Exact){
        Individual_Diagnoses <- Diagnoses_Of_Interest[[Grouping_Name]]
      } else {
        Individual_Diagnoses <- Group %>% group_by(Diagnosis_Name) %>% summarise(.groups = 'drop') %>% pull()
      }
      
      # Look for the individual diagnoses
      for (Diagnosis in Individual_Diagnoses){
        Subgroup_Header <- gsub("_+", "_", gsub(" |\\W", "_", Diagnosis))
        if (Exact){
          Subgroup <- Group %>% filter(Diagnosis_Name == Diagnosis) %>% group_by(EMPI)
        } else {
          Subgroup <- Group %>% filter(grepl(Diagnosis, Diagnosis_Name)) %>% group_by(EMPI)
        }
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
        
        Output_Columns <- Subgroup %>% arrange(Date) %>%
          summarise(!!(as.symbol(Subgroup_Header)) := "Yes",
                    !!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_first_date"))) := first(Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_most_recent_date"))) := last(Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Date, collapse = ";"),
                    .groups = 'drop')
        DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
        DF_to_fill <- DF_to_fill %>% 
          mutate(!!(as.symbol(Subgroup_Header)) := ifelse(is.na(!!(as.symbol(Subgroup_Header))),
                                                          "No", "Yes"))
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
    if(write_files){
      fwrite(Group, str_c(output_file_header, Group_Header, output_file_ending))
    }
    Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
    if (Group_Info){
      Output_Columns <- Group %>% group_by(EMPI) %>% arrange(Date) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_total_diagnoses"))) := n(),
                  !!(as.symbol(str_c(Group_Header, "_first_date"))) := first(Date),
                  !!(as.symbol(str_c(Group_Header, "_most_recent_date"))) := last(Date),
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Date, collapse = ";"),
                  !!(as.symbol(str_c(Group_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                  .groups = 'drop')
      DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
      DF_to_fill <- DF_to_fill %>% select(EMPI:!!(as.symbol(Group_Header)), contains(Group_Header), everything())
    }
    
    if(!missing(Merge_Group_Info_Name)){
      Relevant_rows <- rbind(Relevant_rows, Group)
    }
    rm(Group_Header)
    rm(Group, Output_Columns, Dia_abn)
  }
  if(!missing(Merge_Group_Info_Name)){
    Merge_Header <- str_c("Any_", gsub(" ", "_", Merge_Group_Info_Name))
    Output_Columns <- Relevant_rows %>% arrange(Date) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merge_Header))) := "Yes",
                !!(as.symbol(str_c(Merge_Header, "_total_diagnoses"))) := n(),
                !!(as.symbol(str_c(Merge_Header, "_first_date"))) := first(Date),
                !!(as.symbol(str_c(Merge_Header, "_most_recent_date"))) := last(Date),
                !!(as.symbol(str_c(Merge_Header, "_dates"))) := paste(Date, collapse = ";"),
                !!(as.symbol(str_c(Merge_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                .groups = 'drop')
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
    DF_to_fill <- DF_to_fill %>% mutate(!!(as.symbol(str_c(Merge_Header))) := ifelse(is.na(!!(as.symbol(str_c(Merge_Header)))),
                                                                                     "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Merge_Group_Info_Name), " diagnosis"))
    DF_to_fill <- DF_to_fill %>% select(Original_Columns, starts_with(Merge_Header), everything())
    
    if(write_files){
      fwrite(Relevant_rows, str_c(output_file_header, Merge_Header, output_file_ending))
    }
    rm(Relevant_rows)
  }
  rm(Diagnoses)
  return(DF_to_fill)
}

process_diagnoses_date_cutoff <- function(DF_to_fill_cutoff,
                                          input_file_header_cutoff,
                                          input_file_ending_cutoff,
                                          path_dia_abn_cutoff,
                                          Diagnoses_Of_Interest_cutoff,
                                          Exact_cutoff,
                                          Individual_Info_cutoff,
                                          Group_Info_cutoff,
                                          Merge_Group_Info_Name_cutoff,
                                          write_files_cutoff,
                                          output_file_header_cutoff,
                                          output_file_ending_cutoff,
                                          date_variable_cutoff,
                                          restrict_to_before_cutoff){
  if (missing(Diagnoses_Of_Interest_cutoff)){
    logerror("No list of Diagnoses were specified. Process stopped.")
    return(DF_to_fill_cutoff)
  }
  # Usually this would go in the function but because the compare function doesn't like preassigned information, we do this here
  if (missing(DF_to_fill_cutoff)) {DF_to_fill_cutoff = All_merged}
  if (missing(input_file_header_cutoff)) {input_file_header_cutoff = config$rpdr_file_header}
  if (missing(input_file_ending_cutoff)) {input_file_ending_cutoff = config$rpdr_file_ending}
  if (missing(path_dia_abn_cutoff)) {path_dia_abn_cutoff = str_c(config$data_dir, "Diagnoses_abnormalities/")}
  if (missing(Exact_cutoff)) {Exact_cutoff = FALSE}
  if (missing(Individual_Info_cutoff)) {Individual_Info_cutoff = TRUE}
  if (missing(Group_Info_cutoff)) {Group_Info_cutoff = TRUE}
  if (missing(write_files_cutoff)) {write_files_cutoff = config$create_intermediates}
  if (missing(output_file_header_cutoff)) {output_file_header_cutoff = config$intermediate_files_dir}
  if (missing(output_file_ending_cutoff)) {output_file_ending_cutoff = config$general_file_ending}
  if (missing(restrict_to_before_cutoff)) {restrict_to_before_cutoff = FALSE}
  
  loginfo("Processing diagnoses file...")
  Diagnoses <- data.table(fread(str_c(input_file_header_cutoff, "Dia", input_file_ending_cutoff))) %>% 
    arrange(EMPI, Date)
  if (!dir.exists(path_dia_abn_cutoff)) {dir.create(path_dia_abn_cutoff)}
  path_dia_abn_cutoff <- str_c(path_dia_abn_cutoff, ifelse(restrict_to_before_cutoff,
                                                           "Dia_PreDate_",
                                                           "Dia_PostDate_"))
  # store relevant
  if(!missing(Merge_Group_Info_Name_cutoff)){
    Relevant_rows <- NULL
    Original_Columns <- names(DF_to_fill_cutoff)
  }
  if(write_files_cutoff){
    output_file_header_cutoff <- ifelse(restrict_to_before_cutoff,
                                        str_c(output_file_header_cutoff, "Dia_PreDate_"),
                                        str_c(output_file_header_cutoff, "Dia_PostDate_"))
  }
  loginfo("Restricting prescription by cutoff...")
  # Restrict ids to only the ones in the reduced id list
  Diagnoses <- Diagnoses %>% filter(EMPI %in% DF_to_fill_cutoff$EMPI)
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
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
  # Restrict Diagnoses to only before or only after
  Diagnoses <- left_join(Diagnoses, EMPI_Date_Limit, by = "EMPI")
  time_str = ifelse(restrict_to_before_cutoff, "before", "after")
  loginfo(str_c("Restricting diagnoses to ", time_str, " cutoff dates... "))
  
  Original_row_length <- nrow(Diagnoses)
  if(restrict_to_before_cutoff){
    Diagnoses <- Diagnoses %>% filter(Date <= Cutoff_Date)
  } else {
    Diagnoses <- Diagnoses %>% filter(Date >= Cutoff_Date)
  }
  New_row_length <- nrow(Diagnoses)
  logdebug(str_c(Original_row_length, " diagnoses reduced to ", New_row_length, " diagnoses"))
  rm(time_str, Original_row_length, New_row_length)
  
  # Get the "Any exist" first to lower the search group/increase speed later
  for (Grouping_Name in names(Diagnoses_Of_Interest_cutoff)){
    Group_Header = str_c(ifelse(restrict_to_before_cutoff, "Dia_PreDate_", "Dia_PostDate_"),
                         "Any_", gsub(" ", "_", Grouping_Name))
    Group <- Diagnoses %>%
      filter(grepl(str_c(Diagnoses_Of_Interest_cutoff[[Grouping_Name]], collapse = "|"), Diagnosis_Name))
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Diagnosis_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes",
                .groups = 'drop')
    if (Group_Info_cutoff){
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      DF_to_fill_cutoff <- DF_to_fill_cutoff %>%
        mutate(!!(as.symbol(Group_Header)) := ifelse(is.na(!!(as.symbol(Group_Header))), "No", "Yes"))
    }
    loginfo(str_c(nrow(Output_Columns), " subjects have any ", Grouping_Name, " diagnosis",
                  ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction"))
    rm(Output_Columns)
    
    if (Individual_Info_cutoff){
      if(Exact_cutoff){
        Individual_Diagnoses <- Diagnoses_Of_Interest_cutoff[[Grouping_Name]]
      } else {
        Individual_Diagnoses <- Group %>% group_by(Diagnosis_Name) %>% summarise(.groups = 'drop') %>% pull()
      }
      
      # Look for the individual diagnoses
      for (Diagnosis in Individual_Diagnoses){
        Subgroup_Header <- str_c(ifelse(restrict_to_before_cutoff, "Dia_PreDate_", "Dia_PostDate_"),
                                 gsub("_+", "_", gsub(" |\\W", "_", Diagnosis)))
        if (Exact_cutoff){
          Subgroup <- Group %>% filter(Diagnosis_Name == Diagnosis) %>% group_by(EMPI)
        } else {
          Subgroup <- Group %>% filter(grepl(Diagnosis, Diagnosis_Name)) %>% group_by(EMPI)
        }
        Dia_abn <- Subgroup[duplicated(Subgroup) | duplicated(Subgroup, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " completely duplicated row(s) out of ",
                        nrow(Subgroup), " found. Duplicates removed."))
          fwrite(Dia_abn, str_c(path_dia_abn_cutoff, "Abnormality_1_Duplicate_rows_",
                                Subgroup_Header, output_file_ending_cutoff))
          Subgroup <- Subgroup %>% unique()
        }
        Subgroup2 <- Subgroup %>% select(EMPI, Date, Diagnosis_Name, Hospital)
        Dia_abn <- Subgroup[duplicated(Subgroup2) | duplicated(Subgroup2, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Subgroup),
                        " found. ", sum(duplicated(Subgroup2, fromLast = TRUE)), " rows removed."))
          fwrite(Dia_abn, str_c(path_dia_abn_cutoff, "Abnormality_2_Duplicate_rows_",
                                Subgroup_Header, output_file_ending_cutoff))
          Subgroup <- distinct(Subgroup, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup2)
        
        Output_Columns <- Subgroup %>% arrange(Date) %>%
          summarise(!!(as.symbol(Subgroup_Header)) := "Yes",
                    !!(as.symbol(str_c(Subgroup_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Subgroup_Header, "_first_date"))) := first(Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_most_recent_date"))) := last(Date),
                    !!(as.symbol(str_c(Subgroup_Header, "_dates"))) := paste(Date, collapse = ";"),
                    .groups = 'drop')
        DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
        DF_to_fill_cutoff <- DF_to_fill_cutoff %>% 
          mutate(!!(as.symbol(Subgroup_Header)) := ifelse(is.na(!!(as.symbol(Subgroup_Header))),
                                                          "No", "Yes"))
        loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Diagnosis), " diagnosis",
                      ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction"))
        
        if(write_files_cutoff){
          fwrite(Subgroup, str_c(output_file_header_cutoff, Subgroup_Header, output_file_ending_cutoff))
        }
        rm(Subgroup_Header)
        rm(Subgroup, Output_Columns, Dia_abn)
      }
    }
    if (Group_Info_cutoff){
      if (Individual_Info_cutoff){
        Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
      } else {
        Group <- Group %>% unique()
        Group2 <- Group %>% select(EMPI, Date, Diagnosis_Name, Hospital)
        Dia_abn <- Group[duplicated(Group2) | duplicated(Group2, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Group),
                        " found. ", sum(duplicated(Group2, fromLast = TRUE)), " rows removed."))
          fwrite(Dia_abn, str_c(path_dia_abn_cutoff, "Abnormality_2_Duplicate_rows_",
                                Group_Header, output_file_ending_cutoff))
          Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
        }
        rm(Group2)
      }
      if(write_files_cutoff){
        fwrite(Group, str_c(output_file_header_cutoff, Group_Header, output_file_ending_cutoff))
      }
      Output_Columns <- Group %>% group_by(EMPI) %>% arrange(Date) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_total_diagnoses"))) := n(),
                  !!(as.symbol(str_c(Group_Header, "_first_date"))) := first(Date),
                  !!(as.symbol(str_c(Group_Header, "_most_recent_date"))) := last(Date),
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Date, collapse = ";"),
                  !!(as.symbol(str_c(Group_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                  .groups = 'drop')
      DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
      DF_to_fill_cutoff <- DF_to_fill_cutoff %>% select(EMPI:!!(as.symbol(Group_Header)), contains(Group_Header), everything())
    }
    if(!missing(Merge_Group_Info_Name_cutoff)){
      Relevant_rows <- rbind(Relevant_rows, Group)
    }
    rm(Group_Header)
    rm(Group, Output_Columns, Dia_abn)
  }
  if(!missing(Merge_Group_Info_Name_cutoff)){
    Merge_Header <- str_c(ifelse(restrict_to_before_cutoff, "Dia_PreDate_", "Dia_PostDate_"),
                          "Any_", gsub(" ", "_", Merge_Group_Info_Name_cutoff))
    Output_Columns <- Relevant_rows %>% arrange(Date) %>% group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merge_Header))) := "Yes",
                !!(as.symbol(str_c(Merge_Header, "_total_diagnoses"))) := n(),
                !!(as.symbol(str_c(Merge_Header, "_first_date"))) := first(Date),
                !!(as.symbol(str_c(Merge_Header, "_most_recent_date"))) := last(Date),
                !!(as.symbol(str_c(Merge_Header, "_dates"))) := paste(Date, collapse = ";"),
                !!(as.symbol(str_c(Merge_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                .groups = 'drop')
    DF_to_fill_cutoff <- left_join(DF_to_fill_cutoff, Output_Columns, by = "EMPI")
    DF_to_fill_cutoff <- DF_to_fill_cutoff %>% mutate(!!(as.symbol(str_c(Merge_Header))) := ifelse(is.na(!!(as.symbol(str_c(Merge_Header)))),
                                                                                     "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Merge_Group_Info_Name_cutoff), " diagnosis",
                  ifelse(restrict_to_before_cutoff, " before ", " after "), "cutoff restriction"))
    DF_to_fill_cutoff <- DF_to_fill_cutoff %>% select(Original_Columns, starts_with(Merge_Header), everything())
    
    if(write_files_cutoff){
      fwrite(Relevant_rows, str_c(output_file_header_cutoff, Merge_Header, output_file_ending_cutoff))
    }
    rm(Relevant_rows)
  }
  rm(Diagnoses)
  return(DF_to_fill_cutoff)
}

process_diagnoses_date_compare_cutoff <- function(DF_to_fill = ICS_Tested,
                                                  input_file_header = config$rpdr_file_header,
                                                  input_file_ending = config$rpdr_file_ending,
                                                  path_dia_abn = str_c(config$data_dir, "Diagnoses_abnormalities/"),
                                                  Diagnoses_Of_Interest,
                                                  Exact = FALSE,
                                                  Individual_Info = TRUE,
                                                  Group_Info = TRUE,
                                                  Merge_Group_Info_Name,
                                                  write_files = config$create_intermediates,
                                                  output_file_header = config$intermediate_files_dir,
                                                  output_file_ending = config$general_file_ending,
                                                  cutoff_variable){
  DF_to_fill <- process_diagnoses_date_cutoff(DF_to_fill_cutoff = DF_to_fill,
                                              input_file_header_cutoff = input_file_header,
                                              input_file_ending_cutoff = input_file_ending,
                                              path_dia_abn_cutoff = path_dia_abn,
                                              Diagnoses_Of_Interest_cutoff = Diagnoses_Of_Interest,
                                              Exact_cutoff = Exact,
                                              Individual_Info_cutoff = Individual_Info,
                                              Group_Info_cutoff = Group_Info,
                                              Merge_Group_Info_Name_cutoff = Merge_Group_Info_Name,
                                              write_files_cutoff = write_files,
                                              output_file_header_cutoff = output_file_header,
                                              output_file_ending_cutoff = output_file_ending,
                                              date_variable_cutoff = cutoff_variable,
                                              restrict_to_before_cutoff = TRUE)
  DF_to_fill <- process_diagnoses_date_cutoff(DF_to_fill_cutoff = DF_to_fill,
                                              input_file_header_cutoff = input_file_header,
                                              input_file_ending_cutoff = input_file_ending,
                                              path_dia_abn_cutoff = path_dia_abn,
                                              Diagnoses_Of_Interest_cutoff = Diagnoses_Of_Interest,
                                              Exact_cutoff = Exact,
                                              Individual_Info_cutoff = Individual_Info,
                                              Group_Info_cutoff = Group_Info,
                                              Merge_Group_Info_Name_cutoff = Merge_Group_Info_Name,
                                              write_files_cutoff = write_files,
                                              output_file_header_cutoff = output_file_header,
                                              output_file_ending_cutoff = output_file_ending,
                                              date_variable_cutoff = cutoff_variable,
                                              restrict_to_before_cutoff = FALSE)
  Pre_string <- "Dia_PreDate_"
  Post_string <- "Dia_PostDate_"
  if (!missing(Merge_Group_Info_Name)){
    Merged_Group_Header <- str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", Merge_Group_Info_Name)))
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
             !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) := 
               ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_total_diagnoses"))),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                             "No Change",
                             "Increased"),
                      ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                             "Decreased",
                             ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) <
                                      !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                    "Increased",
                                    ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) >
                                             !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                           "Decreased",
                                           "No Change")))))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                  " subjects were only diagnosed before their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                  " subjects were only diagnosed after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                  " subjects were diagnosed before and after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Increased")),
                  " subjects had an increase in diagnoses after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Decreased")),
                  " subjects had a decrease in diagnoses after their test"))
    loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "No Change")),
                  " subjects had no change in the number of diagnoses after their test"))
  }
  for (Grouping_Name in names(Diagnoses_Of_Interest)){
    Group_Header = str_c("Any_", gsub(" ", "_", Grouping_Name))
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
               !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) := 
                 ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_total_diagnoses"))),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                               "No Change",
                               "Increased"),
                        ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                               "Decreased",
                               ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) <
                                        !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                      "Increased",
                                      ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) >
                                               !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                             "Decreased",
                                             "No Change")))))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                    " subjects were only diagnosed before their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                    " subjects were only diagnosed after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                    " subjects were diagnosed before and after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Increased")),
                    " subjects had an increase in diagnoses after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Decreased")),
                    " subjects had a decrease in diagnoses after their test"))
      loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "No Change")),
                    " subjects had no change in the number of diagnoses after their test"))
    }

    if (Individual_Info){
      if(Exact){
        Individual_Diagnoses <- Diagnoses_Of_Interest[[Grouping_Name]]
      } else {
        Individual_Diagnoses <- Group %>% group_by(Diagnosis_Name) %>% summarise(.groups = 'drop') %>% pull()
      }
      # Look for the individual diagnoses
      for (Diagnosis in Individual_Diagnoses){
        Subgroup_Header <- gsub("_+", "_", gsub(" |\\W", "_", Diagnosis))
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
                 !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) := 
                   ifelse(is.na(!!as.symbol(str_c(Pre_Header, "_total_diagnoses"))),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                                 "No Change",
                                 "Increased"),
                          ifelse(is.na(!!as.symbol(str_c(Post_Header, "_total_diagnoses"))),
                                 "Decreased",
                                 ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) <
                                          !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                        "Increased",
                                        ifelse(!!as.symbol(str_c(Pre_Header, "_total_diagnoses")) >
                                                 !!as.symbol(str_c(Post_Header, "_total_diagnoses")),
                                               "Decreased",
                                               "No Change")))))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),
                      " subjects were only diagnosed before their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),
                      " subjects were only diagnosed after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(Compare_Header) == "Yes - Both ranges")),
                      " subjects were diagnosed before and after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Increased")),
                      " subjects had an increase in diagnoses after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "Decreased")),
                      " subjects had a decrease in diagnoses after their test"))
        loginfo(str_c(nrow(filter(DF_to_fill, !!as.symbol(str_c(Compare_Header, "_total_diagnoses")) == "No Change")),
                      " subjects had no change in the number of diagnoses after their test"))
        rm(Subgroup_Header)
      }
      rm(Diagnosis)
    }
    rm(Group_Header)
  }
  return(DF_to_fill)
}