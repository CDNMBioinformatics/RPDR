require(data.table) # fread, fwrite
require(dplyr) # mutate, filter, select, ...
require(stringr) # str_c
require(tidyverse)
require(logging)
require(readr)
require(zeallot) # %<-% (multiple variable assignment)
require(sqldf)

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
                                       return.flu = TRUE,
                                       return.VitD = FALSE){
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
  if (return.VitD){
    loginfo("Creating Vitamin D table...")
    Group_Stats <- left_join(left_join(Cleaned_DF %>% group_by(Any_VitaminD_Results,
                                                               Low_VitaminD_Results,
                                                               Optimal_VitaminD_Results,
                                                               High_VitaminD_Results) %>%
                                         summarise (Count_All = n(),
                                                    Percent_All = n()/nrow(Cleaned_DF)*100),
                                       Alive %>% group_by(Any_VitaminD_Results,
                                                          Low_VitaminD_Results,
                                                          Optimal_VitaminD_Results,
                                                          High_VitaminD_Results) %>%
                                         summarise (Count_Alive = n(),
                                                    Percent_Alive = n()/nrow(Alive)*100)),
                             Deceased %>% group_by(Any_VitaminD_Results,
                                                   Low_VitaminD_Results,
                                                   Optimal_VitaminD_Results,
                                                   High_VitaminD_Results) %>%
                               summarise (Count_Deceased = n(),
                                          Percent_Deceased = n()/nrow(Deceased)*100))
    fwrite(Group_Stats, str_c(stats_file_input, "Vitamin_D", stats_file_ending))
  }
}