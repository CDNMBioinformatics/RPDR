require(data.table) # fread, fwrite
require(dplyr) # filter, mutate, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)

####################
# Helper Functions #
####################
DiagnosesHelper_CoreFunctionality <- function(Diagnoses_DF,
                                              Main_DF,
                                              Static_Columns_Vector,
                                              Diagnoses_List,
                                              Path_Abnormality,
                                              Write_Output,
                                              File_Prefix,
                                              File_Suffix,
                                              Gen_Group_Info,
                                              Gen_Individual_Info,
                                              Name_Merge_Group,
                                              is_Exact,
                                              is_All_Before_After){
  Header_Prefix = case_when(is_All_Before_After == "All" ~    "",
                            is_All_Before_After == "Before" ~ "Dia_PreDate_",
                            is_All_Before_After == "After" ~  "Dia_PostDate_")
  Log_Suffix = case_when(is_All_Before_After == "All" ~    " diagnoses",
                         is_All_Before_After == "Before" ~ " diagnoses before cutoff restriction",
                         is_All_Before_After == "After" ~  " diagnoses after cutoff restriction")
  
  if (!is.null(Name_Merge_Group)){
    Relevant_rows <- NULL
  }
  #Get the "Any exist" first to lower the search group/increase speed later
  for (Grouping_Name in names(Diagnoses_List)){
    Group_Header = str_c(Header_Prefix, "Any_", gsub(" ", "_", Grouping_Name))
    Group <- Diagnoses_DF %>% filter(grepl(str_c(Diagnoses_List[[Grouping_Name]], collapse = "|"), Diagnosis_Name))
    Output_Columns <- Group %>% group_by(EMPI) %>% select(EMPI, Diagnosis_Name) %>% unique() %>%
      summarise(!!(as.symbol(Group_Header)) := "Yes", .groups = 'drop')
    # Doing this first makes column sorting easier
    if (Gen_Group_Info){
      Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")
      Main_DF <- Main_DF %>% mutate(!!as.symbol(Group_Header) := ifelse(is.na(!!as.symbol(Group_Header)), "No", "Yes"))
    }
    loginfo(str_c(nrow(Output_Columns), " subjects have any ", Grouping_Name, Log_Suffix))
    rm(Output_Columns)
    
    if (Gen_Individual_Info){
      if (is_Exact){ # Use explicit list
        Individual_Diagnoses <- Diagnoses_List[[Grouping_Name]]
      } else { # Use list based on regex strings
        Individual_Diagnoses <- Group %>% group_by(Diagnosis_Name) %>% summarise(.groups = 'drop') %>% pull()
      }
      for (Diagnosis in Individual_Diagnoses){
        Individual_Header <- str_c(Header_Prefix, gsub("_+", "_", gsub(" |\\W", "_", Diagnosis)))
        
        # Look for the individual diagnoses
        if (is_Exact){
          Individual <- Diagnoses_DF %>% filter(Diagnosis_Name == Diagnosis) %>% group_by(EMPI)
        } else {
          Individual <- Diagnoses_DF %>% filter(grepl(Diagnosis, Diagnosis_Name)) %>% group_by(EMPI)
        }
        # Full duplicates check
        Dia_abn <- Individual[duplicated(Individual) | duplicated(Individual, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " completely duplicated row(s) out of ",
                        nrow(Individual), " found. Duplicates removed."))
          fwrite(Dia_abn, str_c(Path_Abnormality, "Abnormality_1_Duplicate_rows_",
                                Individual_Header, File_Suffix))
          Individual <- Individual %>% unique()
        }
        rm(Dia_abn)
        # Partial duplicates check (Multiple instances in the same day logged for clerical reasons)
        Subgroup <- Individual %>% select(EMPI, Date, Diagnosis_Name, Hospital)
        Dia_abn <- Individual[duplicated(Subgroup) | duplicated(Subgroup, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Individual),
                        " found. ", sum(duplicated(Subgroup, fromLast = TRUE)), " rows removed."))
          fwrite(Dia_abn, str_c(Path_Abnormality, "Abnormality_2_Duplicate_rows_",
                                Individual_Header, File_Suffix))
          Individual <- distinct(Individual, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup, Dia_abn)
        
        Output_Columns <- Individual %>% arrange(Date) %>%
          summarise(!!(as.symbol(Individual_Header)) := "Yes",
                    !!(as.symbol(str_c(Individual_Header, "_total_dates"))) := n(),
                    !!(as.symbol(str_c(Individual_Header, "_first_date"))) := first(Date),
                    !!(as.symbol(str_c(Individual_Header, "_most_recent_date"))) := last(Date),
                    !!(as.symbol(str_c(Individual_Header, "_dates"))) := paste(Date, collapse = ";"),
                    .groups = 'drop')
        Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")
        Main_DF <- Main_DF %>% 
          mutate(!!(as.symbol(Individual_Header)) := ifelse(is.na(!!(as.symbol(Individual_Header))), "No", "Yes"))
        loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Diagnosis), Log_Suffix))
        
        if(Write_Output){
          fwrite(Individual, str_c(File_Prefix, Individual_Header, File_Suffix))
        }
        rm(Individual, Output_Columns, Individual_Header)
      }
      rm(Individual_Diagnoses)
    }
    if (Gen_Group_Info){
      if (Gen_Individual_Info){ # Already noted abnormalities
        Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
      } else { # Individual abnormalities haven't been recorded so they should be noted in an overall group
        Dia_abn <- Group[duplicated(Group) | duplicated(Group, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " completely duplicated row(s) out of ",
                        nrow(Group), " found. Duplicates removed."))
          fwrite(Dia_abn, str_c(Path_Abnormality, "Abnormality_1_Duplicate_rows_",
                                Group_Header, File_Suffix))
          Group <- Group %>% unique()
        }
        rm(Dia_abn)
        # Partial duplicates check (Multiple instances in the same day logged for clerical reasons)
        Subgroup <- Group %>% select(EMPI, Date, Diagnosis_Name, Hospital)
        Dia_abn <- Group[duplicated(Subgroup) | duplicated(Subgroup, fromLast=TRUE),]
        if (nrow(Dia_abn) > 0){
          logwarn(str_c(nrow(Dia_abn), " partially duplicated row(s) out of ", nrow(Group),
                        " found. ", sum(duplicated(Subgroup, fromLast = TRUE)), " rows removed."))
          fwrite(Dia_abn, str_c(Path_Abnormality, "Abnormality_2_Duplicate_rows_",
                                Group_Header, File_Suffix))
          Group <- distinct(Group, EMPI, Date, Diagnosis_Name, Hospital, .keep_all = TRUE)
        }
        rm(Subgroup, Dia_abn)
      }
      if (Write_Output) {
        fwrite(Group, str_c(File_Prefix, Group_Header, File_Suffix))
      }
      Output_Columns <- Group %>% group_by(EMPI) %>% arrange(Date) %>%
        summarise(!!(as.symbol(str_c(Group_Header, "_total_diagnoses"))) := n(),
                  !!(as.symbol(str_c(Group_Header, "_first_date"))) := first(Date),
                  !!(as.symbol(str_c(Group_Header, "_most_recent_date"))) := last(Date),
                  !!(as.symbol(str_c(Group_Header, "_dates"))) := paste(Date, collapse = ";"),
                  !!(as.symbol(str_c(Group_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                  .groups = 'drop')
      Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")
      Main_DF <- Main_DF %>% select(EMPI:!!(as.symbol(Group_Header)), contains(Group_Header), everything())
      rm(Output_Columns)
    }
    if (!is.null(Name_Merge_Group)){
      Relevant_rows <- rbind(Relevant_rows, Group)
    }
    rm(Group_Header, Group)
  }
  if (!is.null(Name_Merge_Group)){
    Merge_Header <- str_c(Header_Prefix, "Any_", gsub(" ", "_", Name_Merge_Group))
    Output_Columns <- Diagnoses_DF %>%
      arrange(Date) %>%
      group_by(EMPI) %>%
      summarise(!!(as.symbol(str_c(Merge_Header))) := "Yes",
                !!(as.symbol(str_c(Merge_Header, "_total_diagnoses"))) := n(),
                !!(as.symbol(str_c(Merge_Header, "_first_date"))) := first(Date),
                !!(as.symbol(str_c(Merge_Header, "_most_recent_date"))) := last(Date),
                !!(as.symbol(str_c(Merge_Header, "_dates"))) := paste(Date, collapse = ";"),
                !!(as.symbol(str_c(Merge_Header, "_specific_diagnoses"))) := paste(Diagnosis_Name, collapse = ";"),
                .groups = 'drop')
    Main_DF <- left_join(Main_DF, Output_Columns, by = "EMPI")
    Main_DF <- Main_DF %>%
      mutate(!!(as.symbol(str_c(Merge_Header))) := ifelse(is.na(!!(as.symbol(str_c(Merge_Header)))), "No", "Yes"))
    loginfo(str_c(nrow(Output_Columns), " subjects have a(n) ", tolower(Name_Merge_Group), Log_Suffix))
    Main_DF <- Main_DF %>% select(Static_Columns_Vector, starts_with(Merge_Header), everything())
    
    if(Write_Output){
      fwrite(Diagnoses_DF, str_c(File_Prefix, Merge_Header, File_Suffix))
    }
    return(Main_DF)
  }
}

DiagnosesHelper_CompareOutput <- function(Header,
                                          Main_DF){
  Compare_Header <- str_c("Compare_", Header)
  Pre_Header <- str_c(Pre_string, Header)
  Post_Header <- str_c(Post_string, Header)
  Compare_Header_TD <- str_c(Compare_Header, "_total_diagnoses")
  Pre_Header_TD <- str_c(Pre_Header, "_total_diagnoses")
  Post_Header_TD <- str_c(Post_Header, "_total_diagnoses")
  Main_DF <- Main_DF %>%	
    mutate(!!as.symbol(Compare_Header) := 
             case_when(!!as.symbol(Pre_Header) == "Yes" & !!as.symbol(Post_Header) == "Yes" ~ "Yes - Both ranges",
                       !!as.symbol(Pre_Header) == "Yes" & !!as.symbol(Post_Header) != "Yes" ~ "Yes - Predates cutoff",
                       !!as.symbol(Pre_Header) != "Yes" & !!as.symbol(Post_Header) == "Yes" ~ "Yes - Postdates cutoff",
                       !!as.symbol(Pre_Header) != "Yes" & !!as.symbol(Post_Header) != "Yes" ~ "No"),
           !!as.symbol(Compare_Header_TD) :=
             case_when(
               # No info available
               is.na(!!as.symbol(Pre_Header_TD)) & is.na(!!as.symbol(Post_Header_TD)) ~ "No Change",
               # Info only after cutoff or only before cutoff
               is.na(!!as.symbol(Pre_Header_TD)) & !is.na(!!as.symbol(Post_Header_TD)) ~ "Increased",
               !is.na(!!as.symbol(Pre_Header_TD)) & is.na(!!as.symbol(Post_Header_TD)) ~ "Decreased",
               # Info available on both before and after cutoff
               !is.na(!!as.symbol(Pre_Header_TD)) & !is.na(!!as.symbol(Post_Header_TD)) &
                 !!as.symbol(Pre_Header_TD) < !!as.symbol(Post_Header_TD) ~ "Increased",
               !is.na(!!as.symbol(Pre_Header_TD)) & !is.na(!!as.symbol(Post_Header_TD)) &
                 !!as.symbol(Pre_Header_TD) > !!as.symbol(Post_Header_TD) ~ "Decreased",
               !is.na(!!as.symbol(Pre_Header_TD)) & !is.na(!!as.symbol(Post_Header_TD)) &
                 !!as.symbol(Pre_Header_TD) == !!as.symbol(Post_Header_TD) ~ "No Change"))
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header) == "Yes - Predates cutoff")),	
                " subjects were only diagnosed before their test"))	
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header) == "Yes - Postdates cutoff")),	
                " subjects were only diagnosed after their test"))	
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header) == "Yes - Both ranges")),	
                " subjects were diagnosed before and after their test"))	
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header_TD) == "Increased")),	
                " subjects had an increase in diagnoses after their test"))	
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header_TD) == "Decreased")),	
                " subjects had a decrease in diagnoses after their test"))	
  loginfo(str_c(nrow(filter(Main_DF, !!as.symbol(Compare_Header_TD) == "No Change")),	
                " subjects had no change in the number of diagnoses after their test"))
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
#' @param path_dia_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Diagnoses_abnormalities/)
#' @param Diagnoses_of_interest A list of Diagnoses to select to be provided. Each Key is the Group (usually biobank folder) and Values are Individuals (usually subfolders or items).
#' @param Exact Are you using exact search terms or regular expressions (Default = FALSE)
#' @param Individual_Info Do you want information on the Values? (Default = TRUE)
#' @param Group_Info Do you want information on the Keys? (Default = TRUE)
#' @param Merge_Group_Info_Name If included, all Keys merged together into one group
#' @param write_files Do you want to write intermediates to output? (Default = config$create_intermediates)
#' @param output_file_header The path to intermediates outputs (Default = config$intermediate_files_dir)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_diagnoses(Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("Corticoadrenal insufficiency", "Primary adrenocortical insufficiency", "Other adrenocortical insufficiency", "Unspecified adrenocortical insufficiency")))
#' process_diagnoses(Individual_Info = FALSE, Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("adren.* insufficiency$")))
#' process_diagnoses(Merge_Group_Info_Name = Asthma, Individual_Info = FALSE, Group_Info = FALSE, Diagnoses_Of_Interest = list("Mild Asthma" = c("Mild Asthma"), "Severe Asthma" = c("Severe Asthma")))
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
  loginfo("Processing Diagnoses...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(Diagnoses_Of_Interest)) {
    logerror("No list of Diagnoses were specified. Process stopped.")
    return(DF_to_fill)
  }
  Diagnosis_file_name <- str_c(input_file_header, "Dia", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Diagnosis_file_name))
  
  loginfo("List of Diagnoses of interest (formatted Group: Individual_1, ..., Individual_n")
  for (DN in names(Diagnoses_Of_Interest)){
    loginfo(str_c(DN, ": ", str_c(Diagnoses_Of_Interest[[DN]], collapse = ",")))
  }
  rm(DN)
  
  loginfo(str_c("Generate columns for categorization based on multiple different diagnoses (i.e. multiple biobank folders): ", !missing(Merge_Group_Info_Name)))
  if (!missing(Merge_Group_Info_Name)) {
    loginfo(str_c("Name of collective group: ", Merge_Group_Info_Name))
  } else {
    Merge_Group_Info_Name = NULL
  }
  loginfo(str_c("Generate columns for categorization based on individual general diagnoses (i.e. biobank folder): ", Group_Info))
  loginfo(str_c("Generate columns for categorization based on individual specific diagnoses (i.e. biobank subfolder/item): ", Individual_Info))
  if (Individual_Info) {
    loginfo("Use exact terms or regular expressions to find specific diagnoses: ", Exact)
  }
  
  if (!dir.exists(path_dia_abn)) {dir.create(path_dia_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_dia_abn))
  path_dia_abn <- str_c(path_dia_abn, "Dia_")
  loginfo(str_c("Format of abnormality files: ", path_dia_abn, "{abnormality_description}_{category}", output_file_ending))
  
  loginfo(str_c("Write intermediate files: ", write_files))
  if (write_files){
    output_file_header <- str_c(output_file_header, "Dia_")
    loginfo(str_c("Format of intermediate Files: ", output_file_header, "{intermediate_file_description}", output_file_ending))
  }
  
  #########################
  # Actual Function Start #
  #########################
  Original_Columns <- names(DF_to_fill)
  loginfo("Reading file...")
  Diagnoses <- fread(Diagnosis_file_name)
  loginfo("Read complete")
  rm(Diagnosis_file_name)
  # Should already be formatted that way but adjust in case RPDR changes formatting on their end
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  DF_to_fill <- DiagnosesHelper_CoreFunctionality(
    Diagnoses_DF = Diagnoses,
    Main_DF = DF_to_fill,
    Static_Columns_Vector = Original_Columns,
    Diagnoses_List = Diagnoses_Of_Interest,
    Path_Abnormality = path_dia_abn,
    Write_Output = write_files,
    File_Prefix = output_file_header,
    File_Suffix = output_file_ending,
    Gen_Group_Info = Group_Info,
    Gen_Individual_Info = Individual_Info,
    Name_Merge_Group = Merge_Group_Info_Name,
    is_Exact = Exact,
    is_All_Before_After = "All")
  loginfo("Processing complete")
  return(DF_to_fill)
}

#' Read in file with diagnosis information and add relevant columns to existing data frame
#' Information is limited to all before or all after
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_dia_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Diagnoses_abnormalities/)
#' @param Diagnoses_of_interest A list of Diagnoses to select to be provided. Each Key is the Group (usually biobank folder) and Values are Individuals (usually subfolders or items).
#' @param Exact Are you using exact search terms or regular expressions (Default = FALSE)
#' @param Individual_Info Do you want information on the Values? (Default = TRUE)
#' @param Group_Info Do you want information on the Keys? (Default = TRUE)
#' @param Merge_Group_Info_Name If included, all Keys merged together into one group
#' @param write_files Do you want to write intermediates to output? (Default = config$create_intermediates)
#' @param output_file_header The path to intermediates outputs (Default = config$intermediate_files_dir)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param cutoff_variable The character value representing an existing column with date information
#' @param restrict_to_before_cutoff Do you want to restrict information to before [min - cutoff; TRUE] or after [cutoff - max; FALSE] (Default = TRUE)
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_diagnoses_date_cutoff(Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("Corticoadrenal insufficiency", "Primary adrenocortical insufficiency", "Other adrenocortical insufficiency", "Unspecified adrenocortical insufficiency")), cutoff_variable = "Plasma_Date")
#' process_diagnoses_date_cutoff(Individual_Info = FALSE, Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("adren.* insufficiency$")), cutoff_variable = "Plasma_Date")
#' process_diagnoses_date_cutoff(Merge_Group_Info_Name = Asthma, Individual_Info = FALSE, Group_Info = FALSE, Diagnoses_Of_Interest = list("Mild Asthma" = c("Mild Asthma"), "Severe Asthma" = c("Severe Asthma")), cutoff_variable = "Plasma_Date")
process_diagnoses_date_cutoff <- function(DF_to_fill = All_merged,
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
                                          cutoff_variable,
                                          restrict_to_before_cutoff = TRUE){
  loginfo("Processing Diagnoses...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(Diagnoses_Of_Interest)) {
    logerror("No list of diagnoses were specified. Process stopped.")
    return(DF_to_fill)
  }
  if (missing(cutoff_variable)) {
    logerror("No variable was provided for cutoff. Process stopped.")
    return(DF_to_fill)
  }
  Diagnosis_file_name <- str_c(input_file_header, "Dia", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Diagnosis_file_name))
  
  loginfo(str_c("Variable used for date cutoff: ", cutoff_variable))
  loginfo(str_c("Restriction before cutoff or after cutoff: ", ifelse(restrict_to_before_cutoff, "before", "after")))
  cutoff_prefix <- ifelse(restrict_to_before_cutoff, "Dia_PreDate_", "Dia_PostDate_")
  
  loginfo("List of Diagnoses of interest (formatted Group: Individual_1, ..., Individual_n")
  for (DN in names(Diagnoses_Of_Interest)){
    loginfo(str_c(DN, ": ", str_c(Diagnoses_Of_Interest[[DN]], collapse = ",")))
  }
  rm(DN)
  
  loginfo(str_c("Generate columns for categorization based on multiple different diagnoses (i.e. multiple biobank folders): ", !missing(Merge_Group_Info_Name)))
  if (!missing(Merge_Group_Info_Name)) {
    loginfo(str_c("Name of collective group: ", Merge_Group_Info_Name))
  } else {
    Merge_Group_Info_Name = NULL
  }
  loginfo(str_c("Generate columns for categorization based on individual general diagnoses (i.e. biobank folder): ", Group_Info))
  loginfo(str_c("Generate columns for categorization based on individual specific diagnoses (i.e. biobank subfolder/item): ", Individual_Info))
  if (Individual_Info) {
    loginfo("Use exact terms or regular expressions to find specific diagnoses: ", Exact)
  }
  
  if (!dir.exists(path_dia_abn)) {dir.create(path_dia_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_dia_abn))
  path_dia_abn <- str_c(path_dia_abn, cutoff_prefix)
  loginfo(str_c("Format of abnormality files: ", path_dia_abn, "{abnormality_description}_{category}", output_file_ending))
  
  loginfo(str_c("Write intermediate files: ", write_files))
  if (write_files){
    output_file_header <- str_c(output_file_header, cutoff_prefix)
    loginfo(str_c("Format of intermediate Files: ", output_file_header, "{intermediate_file_description}", output_file_ending))
  }
  
  #########################
  # Actual Function Start #
  #########################
  Original_Columns <- names(DF_to_fill)
  loginfo("Reading file...")
  Diagnoses <- fread(Diagnosis_file_name)
  loginfo("Read complete")
  rm(Diagnosis_file_name)
  # Should already be formatted that way but adjust in case RPDR changes formatting on their end
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  # Restrict to ids only in the reduced id list
  Diagnoses <- Diagnoses %>% filter(EMPI %in% DF_to_fill$EMPI)
  # If date variable includes time, cut off the time
  if (!sum(str_count(DF_to_fill %>% pull(cutoff_variable)) == 10) == nrow(DF_to_fill)){
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, cutoff_variable) %>%
      rename(Cutoff_Date_Time = cutoff_variable) %>%
      extract(Cutoff_Date_Time, c("Cutoff_Date", "Cutoff_Time"),
              regex = "(\\d{4}-\\d{2}-\\d{2}) (\\d{2}:\\d{2})", remove = FALSE) %>%
      mutate(Cutoff_Date = ifelse(is.na(Cutoff_Time), Cutoff_Date_Time, Cutoff_Date)) %>%
      select(EMPI, Cutoff_Date)
  } else {
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, cutoff_variable) %>%
      rename(Cutoff_Date = cutoff_variable)
  }
  # Restrict Diagnoses to only before or only after
  Diagnoses <- left_join(Diagnoses, EMPI_Date_Limit, by = "EMPI")
  rm(EMPI_Date_Limit)
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
  DF_to_fill <- DiagnosesHelper_CoreFunctionality(
    Diagnoses_DF = Diagnoses,
    Main_DF = DF_to_fill,
    Static_Columns_Vector = Original_Columns,
    Diagnoses_List = Diagnoses_Of_Interest,
    Path_Abnormality = path_dia_abn,
    Write_Output = write_files,
    File_Prefix = output_file_header,
    File_Suffix = output_file_ending,
    Gen_Group_Info = Group_Info,
    Gen_Individual_Info = Individual_Info,
    Name_Merge_Group = Merge_Group_Info_Name,
    is_Exact = Exact,
    is_All_Before_After = ifelse(restrict_to_before_cutoff, "Before", "After"))
  loginfo("Processing complete")
  return(DF_to_fill)
}

#' Read in file with diagnosis information and add relevant columns to existing data frame
#' Information is limited to a window [min - max]
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_dia_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Diagnoses_abnormalities/)
#' @param Diagnoses_of_interest A list of Diagnoses to select to be provided. Each Key is the Group (usually biobank folder) and Values are Individuals (usually subfolders or items).
#' @param Exact Are you using exact search terms or regular expressions (Default = FALSE)
#' @param Individual_Info Do you want information on the Values? (Default = TRUE)
#' @param Group_Info Do you want information on the Keys? (Default = TRUE)
#' @param Merge_Group_Info_Name If included, all Keys merged together into one group
#' @param write_files Do you want to write intermediates to output? (Default = config$create_intermediates)
#' @param output_file_header The path to intermediates outputs (Default = config$intermediate_files_dir)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param min_dates The character value representing an existing column with date information
#' @param max_dates The character value representing an existing column with date information
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_diagnoses_set_range(Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("Corticoadrenal insufficiency", "Primary adrenocortical insufficiency", "Other adrenocortical insufficiency", "Unspecified adrenocortical insufficiency")), min_dates = "Plasma_Date_First", max_dates = "Plasma_Date_Most_Recent")
#' process_diagnoses_set_range(Individual_Info = FALSE, Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("adren.* insufficiency$")), min_dates = "Plasma_Date_First", max_dates = "Plasma_Date_Most_Recent")
#' process_diagnoses_set_range(Merge_Group_Info_Name = Asthma, Individual_Info = FALSE, Group_Info = FALSE, Diagnoses_Of_Interest = list("Mild Asthma" = c("Mild Asthma"), "Severe Asthma" = c("Severe Asthma")), min_dates = "Plasma_Date_First", max_dates = "Plasma_Date_Most_Recent")
process_diagnoses_set_range <- function(DF_to_fill = All_merged,
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
                                        min_dates,
                                        max_dates){
  loginfo("Processing Diagnoses...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(Diagnoses_Of_Interest)) {
    logerror("No list of diagnoses were specified. Process stopped.")
    return(DF_to_fill)
  }
  if (missing(min_dates)) {
    logerror("No lower bound date variable provided. Process stopped.")
    return(DF_to_fill)
  }
  if (missing(max_dates)) {
    logerror("No upper bound date variable provided. Process stopped.")
    return(DF_to_fill)
  }
  Diagnosis_file_name <- str_c(input_file_header, "Dia", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Diagnosis_file_name))
  
  loginfo(str_c("Reduce EMPIs in diagnoses file: ", restricted_ids))
  
  loginfo("List of Diagnoses of interest (formatted Group: Individual_1, ..., Individual_n")
  for (DN in names(Diagnoses_Of_Interest)){
    loginfo(str_c(DN, ": ", str_c(Diagnoses_Of_Interest[[DN]], collapse = ",")))
  }
  rm(DN)
  
  loginfo(str_c("Generate columns for categorization based on multiple different diagnoses (i.e. multiple biobank folders): ", !missing(Merge_Group_Info_Name)))
  if (!missing(Merge_Group_Info_Name)) {
    loginfo(str_c("Name of collective group: ", Merge_Group_Info_Name))
  } else {
    Merge_Group_Info_Name = NULL
  }
  loginfo(str_c("Generate columns for categorization based on individual general diagnoses (i.e. biobank folder): ", Group_Info))
  loginfo(str_c("Generate columns for categorization based on individual specific diagnoses (i.e. biobank subfolder/item): ", Individual_Info))
  if (Individual_Info) {
    loginfo("Use exact terms or regular expressions to find specific diagnoses: ", Exact)
  }
  
  if (!dir.exists(path_dia_abn)) {dir.create(path_dia_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_dia_abn))
  path_dia_abn <- str_c(path_dia_abn, "Dia_", min_dates, "_to_", max_dates, "_")
  loginfo(str_c("Format of abnormality files: ", path_dia_abn, "{abnormality_description}_{category}", output_file_ending))
  
  loginfo(str_c("Write intermediate files: ", write_files))
  if (write_files){
    output_file_header <- str_c(output_file_header, "Dia_", min_dates, "_to_", max_dates, "_")
    loginfo(str_c("Format of intermediate Files: ", output_file_header, "{intermediate_file_description}", output_file_ending))
  }
  
  #########################
  # Actual Function Start #
  #########################
  Original_Columns <- names(DF_to_fill)
  loginfo("Reading file...")
  Diagnoses <- fread(Diagnosis_file_name)
  loginfo("Read complete")
  rm(Diagnosis_file_name)
  # Should already be formatted that way but adjust in case RPDR changes formatting on their end
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  # Restrict the EMPIs to ones within the range
  loginfo("Restricting to specified date range...")
  logdebug(str_c("Original row length of Diagnoses file: ", nrow(Diagnoses)))
  Date_Range_DF <- DF_to_fill %>% select(EMPI, min_dates, max_dates)
  Diagnoses <- right_join(Diagnoses, Date_Range_DF, by = "EMPI")
  rm(Date_Range_DF)
  logdebug(str_c("Row length of diagnoses file after join to main EMPI list: ", nrow(Diagnoses)))
  Diagnoses <- Diagnoses %>% filter(Date >= get(min_date)) %>% filter(Date <= get(max_dates))
  logdebug(str_c("Row length after restriction to date window: ", nrow(Diagnoses)))
  DF_to_fill <- DiagnosesHelper_CoreFunctionality(
    Diagnoses_DF = Diagnoses,
    Main_DF = DF_to_fill,
    Static_Columns_Vector = Original_Columns,
    Diagnoses_List = Diagnoses_Of_Interest,
    Path_Abnormality = path_dia_abn,
    Write_Output = write_files,
    File_Prefix = output_file_header,
    File_Suffix = output_file_ending,
    Gen_Group_Info = Group_Info,
    Gen_Individual_Info = Individual_Info,
    Name_Merge_Group = Merge_Group_Info_Name,
    is_Exact = Exact,
    is_All_Before_After = "All")
  loginfo("Processing complete")
  return(DF_to_fill)
}

#' Read in file with diagnosis information and add relevant columns to existing data frame
#' Adds in columns for [Before - cutoff], [cutoff - After], and columns comparing ranges
#' 
#' @param DF_to_fill The data frame you want to add columns to (default = All_merged)
#' @param input_file_header The path + prefix of RPDR file to read (default = config$rpdr_file_header)
#' @param input_file_ending The ending of the RPDR file (default = .txt)
#' @param path_dia_abn The path where you want to record the abnormalities in the EMRs (default = {config$data_dir}/Diagnoses_abnormalities/)
#' @param Diagnoses_of_interest A list of Diagnoses to select to be provided. Each Key is the Group (usually biobank folder) and Values are Individuals (usually subfolders or items).
#' @param Exact Are you using exact search terms or regular expressions (Default = FALSE)
#' @param Individual_Info Do you want information on the Values? (Default = TRUE)
#' @param Group_Info Do you want information on the Keys? (Default = TRUE)
#' @param Merge_Group_Info_Name If included, all Keys merged together into one group
#' @param write_files Do you want to write intermediates to output? (Default = config$create_intermediates)
#' @param output_file_header The path to intermediates outputs (Default = config$intermediate_files_dir)
#' @param output_file_ending The ending of output files (default = config$general_file_ending)
#' @param cutoff_variable The character value representing an existing column with date information
#' 
#' @return \code{DF_to_fill} modified with additional columns
#' 
#' @examples
#' process_diagnoses_date_cutoff(Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("Corticoadrenal insufficiency", "Primary adrenocortical insufficiency", "Other adrenocortical insufficiency", "Unspecified adrenocortical insufficiency")), cutoff_variable = "Plasma_Date")
#' process_diagnoses_date_cutoff(Individual_Info = FALSE, Diagnoses_Of_Interest = list("Adrenal insufficiency" = c("adren.* insufficiency$")), cutoff_variable = "Plasma_Date")
#' process_diagnoses_date_cutoff(Merge_Group_Info_Name = Asthma, Individual_Info = FALSE, Group_Info = FALSE, Diagnoses_Of_Interest = list("Mild Asthma" = c("Mild Asthma"), "Severe Asthma" = c("Severe Asthma")), cutoff_variable = "Plasma_Date")
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
  loginfo("Processing Diagnoses...")
  ###############################
  # Setup and Information Steps #
  ###############################
  if (missing(Diagnoses_Of_Interest)) {
    logerror("No list of diagnoses were specified. Process stopped.")
    return(DF_to_fill)
  }
  if (missing(cutoff_variable)) {
    logerror("No variable was provided for cutoff. Process stopped.")
    return(DF_to_fill)
  }
  Diagnosis_file_name <- str_c(input_file_header, "Dia", input_file_ending)
  loginfo(str_c("Name of file to read in: ", Diagnosis_file_name))
  
  loginfo(str_c("Variable used for date cutoff: ", cutoff_variable))
  
  loginfo("List of Diagnoses of interest (formatted Group: Individual_1, ..., Individual_n")
  for (DN in names(Diagnoses_Of_Interest)){
    loginfo(str_c(DN, ": ", str_c(Diagnoses_Of_Interest[[DN]], collapse = ",")))
  }
  rm(DN)
  
  loginfo(str_c("Generate columns for categorization based on multiple different diagnoses (i.e. multiple biobank folders): ", !missing(Merge_Group_Info_Name)))
  if (!missing(Merge_Group_Info_Name)) {
    loginfo(str_c("Name of collective group: ", Merge_Group_Info_Name))
  }
  loginfo(str_c("Generate columns for categorization based on individual general diagnoses (i.e. biobank folder): ", Group_Info))
  loginfo(str_c("Generate columns for categorization based on individual specific diagnoses (i.e. biobank subfolder/item): ", Individual_Info))
  if (Individual_Info) {
    loginfo("Use exact terms or regular expressions to find specific diagnoses: ", Exact)
  }
  
  if (!dir.exists(path_dia_abn)) {dir.create(path_dia_abn, showWarnings = FALSE)}
  loginfo(str_c("Path to abnormalities: ", path_dia_abn))
  loginfo(str_c("Format of abnormality files: ", path_dia_abn, "{cutoff_marker}_{abnormality_description}_{category}", output_file_ending))
  
  loginfo(str_c("Write intermediate files: ", write_files))
  if (write_files){
    loginfo(str_c("Format of intermediate Files: ", "{cutoff_marker}_{intermediate_file_description}", output_file_ending))
  }
  
  #########################
  # Actual Function Start #
  #########################
  Original_Columns <- names(DF_to_fill)
  loginfo("Reading file...")
  Diagnoses <- fread(Diagnosis_file_name)
  loginfo("Read complete")
  rm(Diagnosis_file_name)
  # Should already be formatted that way but adjust in case RPDR changes formatting on their end
  Diagnoses <- Diagnoses %>% mutate(Date = mdy(Date)) %>% arrange(EMPI, Date)
  # Restrict to ids only in the reduced id list
  Diagnoses <- Diagnoses %>% filter(EMPI %in% DF_to_fill$EMPI)
  # If date variable includes time, cut off the time
  if (!sum(str_count(DF_to_fill %>% pull(cutoff_variable)) == 10) == nrow(DF_to_fill)){
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, cutoff_variable) %>%
      rename(Cutoff_Date_Time = cutoff_variable) %>%
      extract(Cutoff_Date_Time, c("Cutoff_Date", "Cutoff_Time"),
              regex = "(\\d{4}-\\d{2}-\\d{2}) (\\d{2}:\\d{2})", remove = FALSE) %>%
      mutate(Cutoff_Date = ifelse(is.na(Cutoff_Time), Cutoff_Date_Time, Cutoff_Date)) %>%
      select(EMPI, Cutoff_Date)
  } else {
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, cutoff_variable) %>%
      rename(Cutoff_Date = cutoff_variable)
  }
  # Restrict Diagnoses to only before or only after
  Diagnoses <- left_join(Diagnoses, EMPI_Date_Limit, by = "EMPI")
  rm(EMPI_Date_Limit)
  loginfo("Restricting diagnoses to cutoff dates... ")
  
  Original_row_length <- nrow(Diagnoses)
  Diagnoses_Before <- Diagnoses %>% filter(Date <= Cutoff_Date)
  Diagnoses_After <- Diagnoses %>% filter(Date >= Cutoff_Date)
  New_row_length <- nrow(Diagnoses)
  logdebug(str_c(Original_row_length, " diagnoses reduced to ", nrow(Diagnoses_Before), " diagnoses in before cutoff"))
  logdebug(str_c(Original_row_length, " diagnoses reduced to ", nrow(Diagnoses_After), " diagnoses in after cutoff"))
  rm(Original_row_length)
  DF_to_fill <- DiagnosesHelper_CoreFunctionality(
    Diagnoses_DF = Diagnoses_Before,
    Main_DF = DF_to_fill,
    Static_Columns_Vector = Original_Columns,
    Diagnoses_List = Diagnoses_Of_Interest,
    Path_Abnormality = path_dia_abn,
    Write_Output = write_files,
    File_Prefix = output_file_header,
    File_Suffix = output_file_ending,
    Gen_Group_Info = Group_Info,
    Gen_Individual_Info = Individual_Info,
    Name_Merge_Group = Merge_Group_Info_Name,
    is_Exact = Exact,
    is_All_Before_After = "Before")
  Original_Columns <- names(DF_to_fill)
  DF_to_fill <- DiagnosesHelper_CoreFunctionality(
    Diagnoses_DF = Diagnoses_After,
    Main_DF = DF_to_fill,
    Static_Columns_Vector = Original_Columns,
    Diagnoses_List = Diagnoses_Of_Interest,
    Path_Abnormality = path_dia_abn,
    Write_Output = write_files,
    File_Prefix = output_file_header,
    File_Suffix = output_file_ending,
    Gen_Group_Info = Group_Info,
    Gen_Individual_Info = Individual_Info,
    Name_Merge_Group = Merge_Group_Info_Name,
    is_Exact = Exact,
    is_All_Before_After = "After")
  loginfo("Processing complete")
  Pre_string <- "Dia_PreDate_"
  Post_string <- "Dia_PostDate_"
  if (!missing(Merge_Group_Info_Name)){
    Merged_Group_Header <- str_c("Any_", gsub("_{1,}", "_", gsub(" |,|-", "_", Merge_Group_Info_Name)))
    DF_to_fill <- DiagnosesHelper_CompareOutput(Header = Merged_Group_Header,
                                                Main_DF = DF_to_fill)
    rm(Merged_Group_Header)
  }
  for (Grouping_Name in names(Diagnoses_Of_Interest)){
    Group_Header = str_c("Any_", gsub(" ", "_", Grouping_Name))
    if (Group_Info){
      DF_to_fill <- DiagnosesHelper_CompareOutput(Header = Group_Header,
                                                  Main_DF = DF_to_fill)
      rm(Group_Header)
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
        DF_to_fill <- DiagnosesHelper_CompareOutput(Header = Subgroup_Header,
                                                    Main_DF = DF_to_fill)
        rm(Subgroup_Header)
      }
      rm(Diagnosis)
    }
  }
  return(DF_to_fill)
}