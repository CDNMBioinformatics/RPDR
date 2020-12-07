# RPDR

General code and functions for processing/cleaning RPDR data lives here as well as the submodules for specific analyses

## List of related Projects

- How to run analysis that cleaned RPDR data for Priya's ICD and Cortisol study: https://changit.bwh.harvard.edu/resta/RPDR_ACTH

- How to run analysis that cleaned RPDR data for Pegasus study: https://changit.bwh.harvard.edu/resta/RPDR_Pegasus

- How to run analysis that cleaned RPDR data for Alberta's GWAS study: https://changit.bwh.harvard.edu/resta/RPDR_GWAS

- How to run analysis that cleaned RPDR data for Covid Positive Plasma study: https://changit.bwh.harvard.edu/resta/RPDR_Covid_Plasma

- How to run analysis that cleaned RPDR data for Covid Vitamin D study: https://changit.bwh.harvard.edu/resta/RPDR_Covid_VitD

- How to run analysis that cleaned RPDR data for Covid Positive study: https://changit.bwh.harvard.edu/resta/RPDR_Covid_Positive

## Script Descriptions

Start_Cleaning.R
- Reads in Demographics and Biobank RPDR files to create initial dataframe for base of cleaned data
- Can add identifiable parameters Zipcode and Employee [if] Con file was selected in RPDR pull

Add_Biobank_Files.R
- Reads in file pulled from Biobank to add additonal parameters not found in RPDR files
- Includes the following functionality if your Biobank file has additional parameters you don't want to include:
  - PPV.NPV.only: Only include the Currated Disease Population Groups from Biobank
  - Asthma.only: Only include the Asthma columns
  - Asthma.COPD.only: Only include the Asthma and COPD named columns
  - Asthma.Tobacco.only: Only include the Asthma, Tobacco, Smoking, or Smoker columns
- Includes the following functionality if you want to clean the data before merging
  - replace.existence.T.F: Change values from "Yes" and "No" to "1" and "0" respectively for collumns named "Existence Yes/No"
  - clean.list: Change the List All Concepts and List All Values in Biobank from "[a1] [a2] ... [an]" to "a1;a2;...;an"
  - get.data.range: Return the time in years between Date_First to Date_Most_Recent

Clean_Diagnoses_File.R
- Reads in Diagnoses RPDR file and adds cleaned *relevant* data to existing data frame
- Includes functionality:
  - Diagnoses_Of_Interest: REQUIRED list of Diagnoses of interest
    - ex: list(DiaGroup1 = c(Dia1, Dia2, Dia3), DiaGroup2 = c(Dia4))
    - Naming list items is only imporant if Group_Info = TRUE
  - Individial_Info: Return existence of diagnosis, total diagnoses, and ordered dates of diagnoses
  - Group_Info: Group multiple diagnoses together and return information on existence, total diagnoses/date, and ordered dates of diagnoses
    - note multiple different diagnoses are on the same date, only 1 date is listed
  - write_files: Return the subfile of the selected diagnoses from Dia RPDR file

Clean_Encounters_File.R

Clean_Health_History_File.R
- Reads in Health History file and adds cleaned *relevant* data to existing data frame
- Includes functionality to examine BMI information:
  - Return_BMI: TRUE if processing is required
  - Specify parameters *Underweight_Normal*, *Normal_Overweight*, and *Overweight_Obese* to create categories:
    - Underweight: [0 <= x < *Underweight_Normal*]
    - Normal: [*Underweight_Normal* <= x < *Normal_Overweight*]
    - Overweight: [*Normal_Overweight* <= x < *Overweight_Obese*]
    - Obese: [*Overweight_Obese* <= x]
- If Return_BMI is TRUE, returns:
  - Generally: existence, number of values available, average, min, max, and median BMI, list all values by sorted by date, list all dates, and list the predicted category based on median
  - Categorically: existence, number of values recorded, all values sorted by date, and all dates
- Includes functionality to examine if influenza vaccine was administered
  - Return_Influenza: TRUE if processing is required
- If Return_Influenza is TRUE, returns:
  - existence, number of vaccinations recorded, most recent date of vaccination, and sorted dates of vaccinations

Clean_Labs_File.R
- process_ACTH_labs: Reads in Labs RPDR file and adds cleaned ACTH, Cortisol, and DHEA(s) laboratory data to existing data frame
- process_IGE_IGG_labs: Reads in Labs RPDR file and adds cleaned IGE and IGG laboratory data to existing data frame
- process_Covid_labs: Reads in Labs RPDR file and adds cleaned Covid19 testing laboratory data to existing data frame
- process_Cholesterol_labs: Reads in Labs RPDR file and adds cleaned cholesterol laboratory data to existing data frame
- process_VitD_labs: Reads in Labs RPDR file and adds cleaned Vitamin D laboratory data to existing data frame
  - Generally and per category:returns existence, total number of results, median, min, max, and mean Vitamin D level, return sorted dates and values by sorted dates
  - If strict variable is TRUE,
    - removes all non numeric values (including of <x and >x values)
    - doesn't update blank reference range if only included in notes column

Clean_Medications_File.R
- Reads in Medications RPDR file and adds cleaned *relevant* data to existing data frame
- Includes functionality:
  - Group_Column: select method for classifying medications based on Medication_Mapping.txt
  - Medications_Of_Interest: REQUIRED list of Medications of interest
    - ex: list(MedGroup1 = c(Med1, Med2, Med3), MedGroup2 = c(Med4))
    - Naming list items is only imporant if Group_Info = TRUE
  - Individial_Info: Return existence of medication, total prescriptions, first and last dates of prescription, all sorted prescription dates, the name and count of the most common specific prescription
  - Group_Info: Group multiple medications together and return information on existence, total prescription dates, first and last dates of prescription, ordered dates of prescriptions, the most common prescription type (ex Med1) and count, the most common specific prescription and count
    - note if multiple different medications are on the same date, only 1 date is listed
    - note it is possible that the most common specifc prescription is not from the most common prescription group (see most common prescription for individual if you want that info from the most common group)
  - merged_group_name: If included, create merged info based on all medication groups
  - Daily_Dose_Info: Can calculate and add daily dose information based on available dosage/frequency information
    - note most dose/frequency information is either unavailable or incomplete
  - write_files: Return the subfile of the selected medications from Med RPDR file
  - nebs: Return information of nebulizers
  
Clean_Procedures_File.R

Create_Summary_Files.R
