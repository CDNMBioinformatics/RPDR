require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)
require(zeallot) # %<-% (multiple variable assignment)

# Helper Functions
Create_ACTH_Cortisol_DHEA_Output_Columns <- function(ACTH_Cortisol_DHEA_Group,
                                                     Group_Header,
                                                     DF_to_fill){
  # Splitting up Output Columns step into substeps because some median date selection doesn't work well with summarise
  OC_pregroup <- ACTH_Cortisol_DHEA_Group %>% group_by(EMPI, Reference_Units, Seq_Date, Seq_Time) %>%
    summarise(nResults_Per_Time = n(),
              MinResult = min(Result),
              MaxResult = max(Result),
              MedianResult = median(Result),
              MeanResult = mean(Result),
              .groups = 'drop') %>%
    arrange(Seq_Time) %>% group_by(EMPI, Reference_Units, Seq_Date) %>%
    summarise(nResults_Per_Date = sum(nResults_Per_Time),
              nTimes_Per_Date = n(),
              FirstMin = first(MinResult),
              FirstMax = first(MaxResult),
              FirstMedian = first(MedianResult),
              FirstMean = first(MeanResult),
              FirstTime = first(Seq_Time),
              .groups = 'drop') %>%
    mutate(Seq_Date_Time = ifelse(is.na(FirstTime), as.character(Seq_Date), str_c(Seq_Date, " ", FirstTime))) %>% 
    group_by(EMPI, Reference_Units) %>%
    rename(!!as.symbol(str_c(Group_Header, "_Reference_Units")) := Reference_Units)
  DF_to_fill <- left_join(DF_to_fill, OC_pregroup %>% summarise(.groups = 'drop'), by = "EMPI")
  OC_pregroup <- OC_pregroup %>% group_by(EMPI) %>% arrange(Seq_Date_Time)
  
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_nTotalDates")) := n(),
              !!as.symbol(str_c(Group_Header, "_nTotalDatesTimes")) := sum(nTimes_Per_Date),
              !!as.symbol(str_c(Group_Header, "_nTotalResults")) := sum(nResults_Per_Date),
              !!as.symbol(str_c(Group_Header, "_All_Seq_Date_Times")) := paste(Seq_Date_Time, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Min_Result")) := min(FirstMin),
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_First")) := Seq_Date_Time[first(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Min_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMin == min(FirstMin)))],
              !!as.symbol(str_c(Group_Header, "_All_Min_Results")) := paste(FirstMin, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Max_Result")) := max(FirstMax),
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_First")) := Seq_Date_Time[first(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Max_Result_Date_Last")) := Seq_Date_Time[last(which(FirstMax == max(FirstMax)))],
              !!as.symbol(str_c(Group_Header, "_All_Max_Results")) := paste(FirstMax, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  Output_Columns <- OC_pregroup %>% filter(abs(median(FirstMedian) - FirstMedian) == median(FirstMedian) - FirstMedian) %>%
    arrange(EMPI, FirstMedian) %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Median_Result")) := max(FirstMedian),
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_First_or_closest_below")) := Seq_Date_Time[first(which(FirstMedian == max(FirstMedian)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Median_Result_Date_Last_or_closest_below")) := Seq_Date_Time[last(which(FirstMedian == max(FirstMedian)))],
              !!as.symbol(str_c(Group_Header, "_All_Median_Results")) := paste(FirstMedian, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  # All Median is done separately because the previous computation dumps rows that we don't want dumped for the paste
  Output_Columns <- OC_pregroup %>%
    arrange(EMPI, FirstMean) %>%
    summarise(!!as.symbol(str_c(Group_Header, "_Overall_Mean_Result")) := mean(FirstMean),
              !!as.symbol(str_c(Group_Header, "_Overall_Mean_Result_Date_First_or_closest_below")) := Seq_Date_Time[first(which(FirstMean == max(FirstMean)))],
              !!as.symbol(str_c(Group_Header, "_Overall_Mean_Result_Date_Last_or_closest_below")) := Seq_Date_Time[last(which(FirstMean == max(FirstMean)))],
              !!as.symbol(str_c(Group_Header, "_All_Mean_Results")) := paste(FirstMean, collapse = ";"),
              .groups = 'drop')
  
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  loginfo(str_c(nrow(Output_Columns), " subjects have a ", Group_Header, " test performed with useful information"))
  rm(OC_pregroup, Output_Columns)
  return(DF_to_fill)
}

# Test specific functions
process_ACTH_labs <- function(DF_to_fill = All_merged,
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
  if(write_files){
    output_file_header <- str_c(output_file_header, "Lab_")
  }
  if (!exists("strict")){
    strict = FALSE
  }
  if (!exists("skip_ACTH")){
    skip_ACTH = FALSE
  }
  if (!exists("create_cortisol_group")){
    create_cortisol_group = FALSE
  }
  # Only care about ACTH, Cortisol, and/or DHEA(s) in this function
  Labs <- Labs %>% filter(grepl("ACTH|Cortisol($| |, P)|CRH|DHEA", Group_Id))
  
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
  #   Remove the outliers that match the first criteria but are actually not valid
  Labs <- Labs %>%
    mutate(Result = gsub("<0((\\.|0)*)$", "0", Result),
           Result = gsub("(.*(\\d|\\.)+).*", "\\1", Result),
           LessThanX = as.numeric(ifelse(grepl("<", Result), gsub("<((\\d|\\.)+)", "\\1", Result), NA)),
           GreaterThanX = as.numeric(ifelse(grepl(">", Result), gsub(">((\\d|\\.)+)", "\\1", Result), NA)),
           Result = as.numeric(Result),
           Result = ifelse(is.na(LessThanX), Result, LessThanX - min(LessThanX, na.rm = TRUE) / 10),
           Result = ifelse(is.na(GreaterThanX), Result, GreaterThanX + 1/9)) %>%
    select(-c(LessThanX, GreaterThanX)) %>% filter(!is.na(Result))
  
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
  
  Group_Id_list <- Labs %>% group_by(Group_Id) %>% summarise(.groups = 'drop') %>% pull(Group_Id)
  unit_values <- c(dl = 1e-1, ml = 1e-3, ug = 1e-6, ng = 1e-9, pg = 1e-12)
  if (skip_ACTH) {Group_Id_list <- grep("^(?!ACTH).*$", Group_Id_list, value = TRUE, perl = TRUE)}
  # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
  cushings_threshold <- c(ug_dl = 1e2, ng_dl = 1e5, ug_ml = 1, ng_ml = 1e3)
  
  if(create_cortisol_group){
    if(skip_ACTH){
      Cortisol_Group_Id_list <- grep("^Cortisol((?!ACTH).)*$", Group_Id_list, value = TRUE, perl = TRUE)
    } else {
      Cortisol_Group_Id_list <- grep("Cortisol", Group_Id_list, value = TRUE)
    }
    Original_Columns <- names(DF_to_fill)
  }
  
  for (Id in Group_Id_list){
    header <- gsub("\\)", "", gsub("_+", "_", gsub("( |\\(|/|,)", "_", Id)))
    if (strict) { header <- str_c(header, "_Strict")}
    Subgroup <- Labs %>% filter(Group_Id == Id) %>% group_by(EMPI)
    logdebug(Id)
    logdebug(Subgroup %>% group_by(Reference_Units) %>% summarise(n = n(), .groups = 'drop'))
    
    # Note any possible duplicates, but don't necessarily remove
    nSubjects <- Subgroup %>% group_by(EMPI) %>% summarise(count = n(), .groups = 'drop') %>% pull(count) %>% length()
    logdebug(str_c("Number of Subjects: ", nSubjects))
    select_EMPIs <- Subgroup %>% arrange(EMPI) %>% group_by(EMPI, Seq_Date_Time) %>%
      summarise(Count = n(), .groups = 'drop') %>% filter(Count > 1) %>% pull(EMPI) %>% unique()
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
                        Subgroup %>% group_by(Reference_Units) %>%
                          summarise(count = n(), .groups = 'drop') %>%
                          filter(count == max(count)) %>% pull(Reference_Units))
    c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
    
    # Figure out if any reference range information is given as next few steps require for cleaning
    ref_ranges_summary <- Subgroup %>% group_by(Reference_Range) %>% summarise(.groups = 'drop') %>% pull()
    if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
      logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Subgroup), " entries"))
      max_range = Subgroup %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
    } else {
      option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
      option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
      max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                       Subgroup %>% filter(Reference_Units == main_unit) %>%
                                         group_by(Reference_Range) %>%
                                         summarise(.groups = 'drop') %>% pull())), na.rm = TRUE)
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
    logdebug(str_c("Number of Subjects: ", Cortisol_group %>% group_by(EMPI) %>%
                     summarise(count = n(), .groups = 'drop') %>% pull(count) %>% length()))
    # Cortisol should all be the same unit (usually ug/dl) but if that is not the case, change them to whichever unit is the most common
    if (Cortisol_group %>% group_by(Reference_Units) %>%
        summarise(.groups = 'drop') %>% pull(Reference_Units) %>% length() > 1){
      logdebug(Cortisol_group %>% group_by(Reference_Units) %>% summarise(n = n(), .groups = 'drop'))
      # Find the main reference
      main_unit <- Cortisol_group %>% group_by(Reference_Units) %>%
        summarise(count = n(), .groups = 'drop') %>%
        filter(count == max(count)) %>% pull(Reference_Units)
      c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
      # Figure out if any reference range information is given as next few steps require for cleaning
      ref_ranges_summary <- Cortisol_group %>% group_by(Reference_Range) %>%
        summarise(.groups = 'drop') %>% pull()
      if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
        logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Cortisol_group), " entries"))
        max_range = Cortisol_group %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
      } else {
        option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
        option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
        max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                         Cortisol_group %>% filter(Reference_Units == main_unit) %>%
                                           group_by(Reference_Range) %>% summarise(.groups = 'drop') %>% pull())), na.rm = TRUE)
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
    DF_to_fill <- DF_to_fill %>% select(Original_Columns, starts_with(header), everything())
    rm(Cortisol_group, header)
  }
  rm(Labs)
  return(DF_to_fill)
}

# Test specific functions
process_ACTH_labs_cutoff <- function(DF_to_fill = All_merged,
                                     input_file_header = config$rpdr_file_header,
                                     input_file_ending = config$rpdr_file_ending,
                                     path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                                     skip_ACTH = config$ACTH_params$skip_ACTH,
                                     strict = config$ACTH_params$strict,
                                     create_cortisol_group = config$ACTH_params$create_cortisol_group,
                                     write_files = config$create_intermediates,
                                     output_file_header = config$intermediate_files_dir,
                                     output_file_ending = config$general_file_ending,
                                     date_variable_cutoff,
                                     restrict_to_before_cutoff = FALSE){
  loginfo("Processing labs file...")
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI, Seq_Date_Time)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  path_lab_abn <- str_c(path_lab_abn, ifelse(restrict_to_before_cutoff,
                                             "Lab_PreDate_",
                                             "Lab_PostDate_"))
  if(write_files){
    output_file_header <- ifelse(restrict_to_before_cutoff,
                                 str_c(output_file_header, "Lab_PreDate_"),
                                 str_c(output_file_header, "Lab_PostDate_"))
  }
  if (!exists("strict")){
    strict = FALSE
  }
  if (!exists("skip_ACTH")){
    skip_ACTH = FALSE
  }
  if (!exists("create_cortisol_group")){
    create_cortisol_group = FALSE
  }
  #Restrict ids to reduced list
  Labs <- Labs %>% filter(EMPI %in% DF_to_fill$EMPI)
  if (!sum(str_count(DF_to_fill %>% pull(date_variable_cutoff)) == 10) == nrow(DF_to_fill)){
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, date_variable_cutoff) %>%
      rename(Cutoff_Date_Time = date_variable_cutoff) %>%
      extract(Cutoff_Date_Time, c("Cutoff_Date", "Cutoff_Time"),
              regex = "(\\d{4}-\\d{2}-\\d{2}) (\\d{2}:\\d{2})", remove = FALSE) %>%
      mutate(Cutoff_Date = ifelse(is.na(Cutoff_Time), Cutoff_Date_Time, Cutoff_Date)) %>% select(EMPI, Cutoff_Date)
  } else {
    EMPI_Date_Limit <- DF_to_fill %>% select(EMPI, date_variable_cutoff) %>%
      rename(Cutoff_Date = date_variable_cutoff)
  }
  # Restrict Labs to only before or only after date
  Labs <- left_join(Labs, EMPI_Date_Limit, by = "EMPI")
  
  # Clean data (6 in the main fn but 1 now to make cutoff possible)
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
  time_str = ifelse(restrict_to_before_cutoff, "before", "after")
  loginfo(str_c("Restricting diagnoses to ", time_str, " cutoff dates... "))
  
  Original_row_length <- nrow(Labs)
  if(restrict_to_before_cutoff){
    Labs <- Labs %>% filter(Seq_Date <= Cutoff_Date)
  } else {
    Labs <- Labs %>% filter(Seq_Date >= Cutoff_Date)
  }
  New_row_length <- nrow(Labs)
  logdebug(str_c(Original_row_length, " lab entries reduced to ", New_row_length, " lab entries"))
  rm(time_str, Original_row_length, New_row_length)
  
  
  # Only care about ACTH, Cortisol, and/or DHEA(s) in this function
  Labs <- Labs %>% filter(grepl("ACTH|Cortisol($| |, P)|CRH|DHEA", Group_Id))
  
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
  #   Remove the outliers that match the first criteria but are actually not valid
  Labs <- Labs %>%
    mutate(Result = gsub("<0((\\.|0)*)$", "0", Result),
           Result = gsub("(.*(\\d|\\.)+).*", "\\1", Result),
           LessThanX = as.numeric(ifelse(grepl("<", Result), gsub("<((\\d|\\.)+)", "\\1", Result), NA)),
           GreaterThanX = as.numeric(ifelse(grepl(">", Result), gsub(">((\\d|\\.)+)", "\\1", Result), NA)),
           Result = as.numeric(Result),
           Result = ifelse(is.na(LessThanX), Result, LessThanX - min(LessThanX, na.rm = TRUE) / 10),
           Result = ifelse(is.na(GreaterThanX), Result, GreaterThanX + 1/9)) %>%
    select(-c(LessThanX, GreaterThanX)) %>% filter(!is.na(Result))
  
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
  
  Group_Id_list <- Labs %>% group_by(Group_Id) %>% summarise(.groups = 'drop') %>% pull(Group_Id)
  unit_values <- c(dl = 1e-1, ml = 1e-3, ug = 1e-6, ng = 1e-9, pg = 1e-12)
  if (skip_ACTH) {Group_Id_list <- grep("^(?!ACTH).*$", Group_Id_list, value = TRUE, perl = TRUE)}
  # According to readings, cortisol levels in the 50-100 ug/dl correspond to cushings syndorm which is the max value we want to filter.
  cushings_threshold <- c(ug_dl = 1e2, ng_dl = 1e5, ug_ml = 1, ng_ml = 1e3)
  
  if(create_cortisol_group){
    if(skip_ACTH){
      Cortisol_Group_Id_list <- grep("^Cortisol((?!ACTH).)*$", Group_Id_list, value = TRUE, perl = TRUE)
    } else {
      Cortisol_Group_Id_list <- grep("Cortisol", Group_Id_list, value = TRUE)
    }
    Original_Columns <- names(DF_to_fill)
  }
  
  for (Id in Group_Id_list){
    header <- gsub("\\)", "", gsub("_+", "_", gsub("( |\\(|/|,)", "_", Id)))
    if (strict) { header <- str_c(header, "_Strict")}
    Subgroup <- Labs %>% filter(Group_Id == Id) %>% group_by(EMPI)
    logdebug(Id)
    logdebug(Subgroup %>% group_by(Reference_Units) %>% summarise(n = n(), .groups = 'drop'))
    
    # Note any possible duplicates, but don't necessarily remove
    nSubjects <- Subgroup %>% group_by(EMPI) %>% summarise(count = n(), .groups = 'drop') %>% pull(count) %>% length()
    logdebug(str_c("Number of Subjects: ", nSubjects))
    select_EMPIs <- Subgroup %>% arrange(EMPI) %>% group_by(EMPI, Seq_Date_Time) %>%
      summarise(Count = n(), .groups = 'drop') %>% filter(Count > 1) %>% pull(EMPI) %>% unique()
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
                        Subgroup %>% group_by(Reference_Units) %>%
                          summarise(count = n(), .groups = 'drop') %>%
                          filter(count == max(count)) %>% pull(Reference_Units))
    c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
    
    # Figure out if any reference range information is given as next few steps require for cleaning
    ref_ranges_summary <- Subgroup %>% group_by(Reference_Range) %>% summarise(.groups = 'drop') %>% pull()
    if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
      logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Subgroup), " entries"))
      max_range = Subgroup %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
    } else {
      option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
      option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
      max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                       Subgroup %>% filter(Reference_Units == main_unit) %>%
                                         group_by(Reference_Range) %>%
                                         summarise(.groups = 'drop') %>% pull())), na.rm = TRUE)
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
    logdebug(str_c("Number of Subjects: ", Cortisol_group %>% group_by(EMPI) %>%
                     summarise(count = n(), .groups = 'drop') %>% pull(count) %>% length()))
    # Cortisol should all be the same unit (usually ug/dl) but if that is not the case, change them to whichever unit is the most common
    if (Cortisol_group %>% group_by(Reference_Units) %>%
        summarise(.groups = 'drop') %>% pull(Reference_Units) %>% length() > 1){
      logdebug(Cortisol_group %>% group_by(Reference_Units) %>% summarise(n = n(), .groups = 'drop'))
      # Find the main reference
      main_unit <- Cortisol_group %>% group_by(Reference_Units) %>%
        summarise(count = n(), .groups = 'drop') %>%
        filter(count == max(count)) %>% pull(Reference_Units)
      c(main_unit_num, main_unit_den) %<-% str_split_fixed(main_unit, "/", n = 2)
      # Figure out if any reference range information is given as next few steps require for cleaning
      ref_ranges_summary <- Cortisol_group %>% group_by(Reference_Range) %>%
        summarise(.groups = 'drop') %>% pull()
      if (length(ref_ranges_summary) == 1 && ref_ranges_summary == ""){
        logwarn(str_c("GroupId ", Id, " does not list reference range information in all ", nrow(Cortisol_group), " entries"))
        max_range = Cortisol_group %>% filter(Abnormal_Flag == "") %>% pull(Result) %>% max()
      } else {
        option1 <- "(.*( |-)((\\d|\\.)*)(\\(*a.+)*$)" # if a-b given, select b (note some mention a.m. after) [select \\3 of 1-5]
        option2 <- "(^<((\\d|\\.)+)$)"                # if  <b given, select b [select \\2 of 1-3]
        max_range <- max(as.numeric(gsub(str_c(option1, option2, sep = "|"), "\\3\\7",
                                         Cortisol_group %>% filter(Reference_Units == main_unit) %>%
                                           group_by(Reference_Range) %>% summarise(.groups = 'drop') %>% pull())), na.rm = TRUE)
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
    DF_to_fill <- DF_to_fill %>% select(Original_Columns, starts_with(header), everything())
    rm(Cortisol_group, header)
  }
  rm(Labs)
  return(DF_to_fill)
}

process_IGE_IGG_labs <- function(DF_to_fill = All_merged,
                                 input_file_header = config$rpdr_file_header,
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
                                                        summarise(Count = n(), .groups = 'drop') %>%
                                                        filter(Count > 1) %>% group_by(EMPI) %>%
                                                        summarise(.groups = 'drop') %>% pull())) %>%
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
      summarise(mean_duplicates = mean(Result),
                .groups = 'drop') %>%
      group_by(EMPI, Reference_Units) %>%
      summarise(!!as.symbol(str_c(ID, "_nTotalValues")) := n(),
                !!as.symbol(str_c(ID, "_Mean")) := mean(mean_duplicates),
                !!as.symbol(str_c(ID, "_SD")) := sd(mean_duplicates),
                !!as.symbol(str_c(ID, "_All_Seq_Date_Times")) := paste(Seq_Date_Time, collapse = ";"),
                !!as.symbol(str_c(ID, "_All_Results")) := paste(mean_duplicates, collapse = ";"),
                .groups = 'drop') %>%
      rename(!!as.symbol(str_c(ID, "_Reference_Units")) := Reference_Units)
    DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  }
  
  rm(Labs)
  return(DF_to_fill)
}

process_Covid_labs <- function(DF_to_fill = All_merged,
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
                                       "Result_Order" = paste(Result_Updated, collapse = ";"),
                                       .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Tested = ifelse(is.na(Tested), "No", Tested),
                                      Number_Of_Total_Tests = ifelse(is.na(Number_Of_Total_Tests),
                                                                     0, Number_Of_Total_Tests))
  Subgroup <- Labs %>% filter(grepl("POSITIVE", Result_Updated))
  Output_Columns <- Subgroup %>% summarise("Positive_Result" = "Yes",
                                           "Number_Of_Total_Positive_Tests" = n(),
                                           "First_Positive_Test_Date" = first(Seq_Date_Time),
                                           "All_Positive_Test_Dates" = paste(Seq_Date_Time, collapse = ";"),
                                           .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Positive_Result = ifelse(is.na(Positive_Result), "No", Positive_Result),
                                      Number_Of_Total_Positive_Tests = ifelse(is.na(Number_Of_Total_Positive_Tests),
                                                                              0, Number_Of_Total_Positive_Tests))
  Subgroup <- Labs %>% filter(grepl("NEGATIVE", Result_Updated))
  Output_Columns <- Subgroup %>% summarise("Negative_Result" = "Yes",
                                           "Number_Of_Total_Negative_Tests" = n(),
                                           "First_Negative_Test_Date" = first(Seq_Date_Time),
                                           "All_Negative_Test_Dates" = paste(Seq_Date_Time, collapse = ";"),
                                           .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns, by = "EMPI")
  DF_to_fill <- DF_to_fill %>% mutate(Negative_Result = ifelse(is.na(Negative_Result), "No", Negative_Result),
                                      Number_Of_Total_Negative_Tests = ifelse(is.na(Number_Of_Total_Negative_Tests),
                                                                              0, Number_Of_Total_Negative_Tests))
  
  rm(Labs, Subgroup)
  return(DF_to_fill)
}

process_Cholesterol_labs <- function(DF_to_fill = All_merged,
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
                                                            "Intermediate", "High")),
              .groups = 'drop')
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
              Cholesterol_optimal_dates = paste(Seq_Date_Time, collapse = ";"),
              .groups = 'drop')
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
              Cholesterol_intermediate_dates = paste(Seq_Date_Time, collapse = ";"),
              .groups = 'drop')
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
              Cholesterol_high_dates = paste(Seq_Date_Time, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>%
    mutate(Cholesterol_high = ifelse(is.na(Cholesterol_high), "No", Cholesterol_high),
           Cholesterol_high_number_recorded = ifelse(is.na(Cholesterol_high_number_recorded),
                                                     0, Cholesterol_high_number_recorded))
  rm(Optimal, Intermediate, High, Labs, Output_Columns)
  return(DF_to_fill)
}

process_VitD_labs <- function(DF_to_fill = All_merged,
                              input_file_header = config$rpdr_file_header,
                              input_file_ending = config$rpdr_file_ending,
                              path_lab_abn = str_c(config$data_dir, "Lab_abnormalities/"),
                              output_file_ending = config$general_file_ending,
                              write_files = config$create_intermediates,
                              strict = FALSE){
  loginfo("Processing Vitamin D labs...")
  Labs <- data.table(fread(str_c(input_file_header, "Lab", input_file_ending))) %>% arrange(EMPI)
  if (!dir.exists(path_lab_abn)) {dir.create(path_lab_abn)}
  logdebug(str_c("Note: All Lab abnormalites can be found at ", path_lab_abn))
  if(write_files){
    output_file_header <- str_c(output_file_header, "Lab_")
  }
  
  Labs <- Labs %>% filter(grepl("[Vv]itamin", Group_Id))
  
  # Cleanup Results
  Labs <- Labs %>% mutate(Result_Update = as.numeric(Result))
  if (strict){
    loginfo("Removing all non-numeric Results...")
  } else {
    loginfo("Adjusting partial non-numeric Results to numeric values...")
    Labs <- Labs %>% mutate(Result_Update = ifelse(grepl(">", Result),
                                                   as.numeric(gsub(".*>\\D*((\\d|\\.)+).*", "\\1", Result)) + 0.5,
                                                   Result_Update),
                            Result_Update = ifelse(grepl("<", Result),
                                                   as.numeric(gsub(".*<\\D*((\\d|\\.)+).*", "\\1", Result)) - 0.5,
                                                   Result_Update))
    loginfo("Removing all remaining non-numeric Results...")
  }
  ## Get rid of invalid results
  Labs <- Labs %>% filter(!is.na(Result_Update))
  
  # Cleanup Range Data for check
  ## Part 1
  Labs <- Labs %>% mutate(Reference_Range_Update = Reference_Range)
  Labs <- Labs %>% mutate(Reference_Range_Update_gt = ifelse(grepl(">", Reference_Range_Update),
                                                      as.numeric(gsub(" *> *((\\d|\\.)+).*$", "\\1", Reference_Range_Update)) + 1,
                                                      NA),
                          Reference_Range_Update = ifelse(!is.na(Reference_Range_Update_gt),
                                                   str_c(Reference_Range_Update_gt, "-1000000"),
                                                   Reference_Range_Update),
                          Reference_Range_Update_lt = ifelse(grepl("<", Reference_Range_Update),
                                                      as.numeric(gsub(" *< *((\\d|\\.)+).*$", "\\1", Reference_Range_Update)) - 1,
                                                      NA),
                          Reference_Range_Update = ifelse(!is.na(Reference_Range_Update_lt),
                                                   str_c("0-",Reference_Range_Update_lt),
                                                   Reference_Range_Update),
                          Reference_Range_Update = ifelse(!grepl("-", Reference_Range_Update),
                                                   "",
                                                   Reference_Range_Update)) %>%
    select(-c(Reference_Range_Update_gt, Reference_Range_Update_lt)) %>%
    extract(Reference_Range_Update, c("Reference_Range_Update_Min", "Reference_Range_Update_Max"),
            regex = "(\\S+) *- *(\\S+)", remove = FALSE) %>%
    mutate(Reference_Range_Update_Min = as.numeric(Reference_Range_Update_Min),
           Reference_Range_Update_Max = as.numeric(Reference_Range_Update_Max),
           Reference_Range_Update = ifelse(!is.na(Reference_Range_Update_Min),
                                           str_c(Reference_Range_Update_Min,
                                                 "-",
                                                 Reference_Range_Update_Max),
                                           Reference_Range_Update))
  Median_Range_Min <- median(Labs$Reference_Range_Update_Min, na.rm = TRUE)
  Median_Range_Max <- median(Labs$Reference_Range_Update_Max, na.rm = TRUE)
  Max_Range_Min <- max(Labs$Reference_Range_Update_Min, na.rm = TRUE)
  
  ## Part 2
  if (strict){
    loginfo("Skipping adjusting Reference_Range...")
  } else {
    loginfo("Deriving Reference_Range from Result_Text...")
    ## Part 2a: Find a - b in text
    Labs <- Labs %>%
      mutate(Reference_Range_Text_1 = ifelse(grepl(".*\\D{2}((\\d|\\.)+ *(-|[Tt][Oo]) *(\\d|\\.)+).*", Result_Text),
                                       gsub(".*\\D{2}((\\d|\\.)+ *(-|[Tt][Oo]) *(\\d|\\.)+).*", "\\1", Result_Text),
                                       "")) %>%
      extract(Reference_Range_Text_1, c("Reference_Range_Text_1_Min","separator", "Reference_Range_Text_1_Max"),
              regex = "(\\S+) *(-|[Tt][Oo]) *(\\S+)", remove = FALSE) %>%
      select(-separator) %>%
      mutate(Reference_Range_Text_1_Min = as.numeric(Reference_Range_Text_1_Min),
             Reference_Range_Text_1_Max = as.numeric(Reference_Range_Text_1_Max)) %>%
      mutate(Reference_Range_Text_1 = ifelse(!is.na(Reference_Range_Text_1_Min),
                                             str_c(Reference_Range_Text_1_Min,
                                                   "-",
                                                   Reference_Range_Text_1_Max),
                                             Reference_Range_Text_1),
             Reference_Range_Text_1 = ifelse(Reference_Range_Text_1_Min > Reference_Range_Text_1_Max,
                                             NA, Reference_Range_Text_1),
             Reference_Range_Text_1 = ifelse(Reference_Range_Text_1_Min > Median_Range_Max,
                                             NA, Reference_Range_Text_1))

    ## Part 2b: Find > a in text
    Labs <- Labs %>%
      mutate(Reference_Range_Text_2 = ifelse(grepl(">|Greater than", Result_Text, ignore.case = TRUE),
                                       as.numeric(gsub(".*(>|Greater than) *((\\d|\\.)+).*", "\\2", Result_Text, ignore.case = TRUE)) + 1,
                                       NA),
             Reference_Range_Text_2 = ifelse(!is.na(Reference_Range_Text_2),
                                       str_c(Reference_Range_Text_2, "-1000000"),
                                       ""))
    
    ## Part 2c: Find < a in text
    Labs <- Labs %>%
      mutate(Reference_Range_Text_3 = ifelse(grepl("<|Less than", Result_Text, ignore.case = TRUE),
                                       as.numeric(gsub(".*(<|Less than) *((\\d|\\.)+).*", "\\2", Result_Text, ignore.case = TRUE)) - 1,
                                       NA),
             Reference_Range_Text_3 = ifelse(!is.na(Reference_Range_Text_3),
                                       str_c("0-", Reference_Range_Text_3),
                                       ""))
    
    ## All the remaining Reference_Range missing values cannot be filled with the available Result_Text
    Labs <- Labs %>% 
      extract(Reference_Range, c("Reference_Range_Min", "Reference_Range_Max"),
              regex = "(\\S+) *- *(\\S+)", remove = FALSE) %>%
      mutate(Reference_Range_Min = as.numeric(Reference_Range_Min),
             Reference_Range_Max = as.numeric(Reference_Range_Max))
  }
  # Cleanup Units
  Labs <- Labs %>% mutate(Reference_Units = tolower(Reference_Units))
  
  # Cleanup Abnormality Flags
  loginfo("Adjusting Abnormality Flags...")
  Labs <- Labs %>% mutate(Abnormal_Flag = ifelse(Reference_Range != "",
                                                 ifelse(Result_Update < Reference_Range_Min,
                                                        "L",
                                                        ifelse(Result_Update > Reference_Range_Max,
                                                               "H",
                                                               "")),
                                                 ifelse(Result_Update < Median_Range_Min,
                                                        "L",
                                                        ifelse(Result_Update > Median_Range_Max,
                                                               "H",
                                                               ""))))
  
  # Cleanup Seq_Date_Time
  ## Current plan: Drop time from object and just keep date as unclear if time is important/
  ##               time is not available for all rows
  Labs <- Labs %>%
    extract(Seq_Date_Time, c("Seq_Date", "Seq_Time"),
            regex = "(\\d{2}/\\d{2}/\\d{4}) (\\d{2}:\\d{2})", remove = FALSE) %>%
    mutate(Seq_Date = ifelse(is.na(Seq_Time), Seq_Date_Time, Seq_Date),
           Seq_Date = as.Date(Seq_Date, "%m/%d/%Y")) %>%
    arrange(Seq_Date)
  
  # Generate Outputs
  Output_Columns <- Labs %>% group_by(EMPI) %>% arrange(Seq_Date) %>%
    summarise(Any_VitaminD_Results = "Yes",
              Any_Total_VitaminD_Results = n(),
              Any_VitaminD_Median_Result = median(Result_Update, na.rm = TRUE),
              Any_VitaminD_Min_Result = min(Result_Update, na.rm = TRUE),
              Any_VitaminD_Max_Result = max(Result_Update, na.rm = TRUE),
              Any_VitaminD_Mean_Result = mean(Result_Update, na.rm = TRUE),
              Any_VitaminD_All_Dates = paste(Seq_Date, collapse = ";"),
              Any_VitaminD_All_Results_By_Date = paste(Result, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>% mutate(Any_VitaminD_Results = ifelse(is.na(Any_VitaminD_Results),
                                                                    "No",
                                                                    Any_VitaminD_Results))
  loginfo(str_c(nrow(Output_Columns), " subjects have any Vitamin D information available"))
  
  Output_Columns <- Labs %>% filter(Abnormal_Flag == "L") %>% group_by(EMPI) %>% arrange(Seq_Date) %>%
    summarise(Low_VitaminD_Results = "Yes",
              Low_Total_VitaminD_Results = n(),
              Low_VitaminD_Median_Result = median(Result_Update, na.rm = TRUE),
              Low_VitaminD_Min_Result = min(Result_Update, na.rm = TRUE),
              Low_VitaminD_Max_Result = max(Result_Update, na.rm = TRUE),
              Low_VitaminD_Mean_Result = mean(Result_Update, na.rm = TRUE),
              Low_VitaminD_All_Dates = paste(Seq_Date, collapse = ";"),
              Low_VitaminD_All_Results_By_Date = paste(Result, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>% mutate(Low_VitaminD_Results = ifelse(is.na(Low_VitaminD_Results),
                                                                    "No",
                                                                    Low_VitaminD_Results))
  loginfo(str_c(nrow(Output_Columns), " subjects have low Vitamin D information available"))
  
  Output_Columns <- Labs %>% filter(Abnormal_Flag == "") %>% group_by(EMPI) %>% arrange(Seq_Date) %>%
    summarise(Optimal_VitaminD_Results = "Yes",
              Optimal_Total_VitaminD_Results = n(),
              Optimal_VitaminD_Median_Result = median(Result_Update, na.rm = TRUE),
              Optimal_VitaminD_Min_Result = min(Result_Update, na.rm = TRUE),
              Optimal_VitaminD_Max_Result = max(Result_Update, na.rm = TRUE),
              Optimal_VitaminD_Mean_Result = mean(Result_Update, na.rm = TRUE),
              Optimal_VitaminD_All_Dates = paste(Seq_Date, collapse = ";"),
              Optimal_VitaminD_All_Results_By_Date = paste(Result, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>% mutate(Optimal_VitaminD_Results = ifelse(is.na(Optimal_VitaminD_Results),
                                                                        "No",
                                                                        Optimal_VitaminD_Results))
  loginfo(str_c(nrow(Output_Columns), " subjects have optimal Vitamin D information available"))
  
  Output_Columns <- Labs %>% filter(Abnormal_Flag == "H") %>% group_by(EMPI) %>% arrange(Seq_Date) %>%
    summarise(High_VitaminD_Results = "Yes",
              High_Total_VitaminD_Results = n(),
              High_VitaminD_Median_Result = median(Result_Update, na.rm = TRUE),
              High_VitaminD_Min_Result = min(Result_Update, na.rm = TRUE),
              High_VitaminD_Max_Result = max(Result_Update, na.rm = TRUE),
              High_VitaminD_Mean_Result = mean(Result_Update, na.rm = TRUE),
              High_VitaminD_All_Dates = paste(Seq_Date, collapse = ";"),
              High_VitaminD_All_Results_By_Date = paste(Result, collapse = ";"),
              .groups = 'drop')
  DF_to_fill <- left_join(DF_to_fill, Output_Columns)
  DF_to_fill <- DF_to_fill %>% mutate(High_VitaminD_Results = ifelse(is.na(High_VitaminD_Results),
                                                                     "No",
                                                                     High_VitaminD_Results))
  loginfo(str_c(nrow(Output_Columns), " subjects have high Vitamin D information available"))
  return(DF_to_fill)
}