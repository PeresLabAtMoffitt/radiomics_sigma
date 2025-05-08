# Import library
library(tidyverse)
library(lubridate)


#################################################################### I ### Load data ----
path <- fs::path("","Volumes","Peres_Research", "Ovarian - Radiomics")
ovary_data <-
  read_csv(paste0(path, "/SIGMA",
                  "/Data_05082025/CANCER_REGISTRY_SEROUS_OVARY.csv")) %>% 
  janitor::clean_names()


#################################################################### II ### Data cleaning ----
# Filter for HGSC
table(ovary_data$cdsc_primary_site_group, useNA = "always")
table(ovary_data$primary_site_cd, useNA = "always")
table(ovary_data$histology_desc, useNA = "always")

# Check treatment variables
which(ovary_data$summary_of_rx_1st_course_desc != ovary_data$summary_of_rx_1st_course_at_this_hospital_desc)
ovary_data$summary_of_rx_1st_course_desc[ovary_data$summary_of_rx_1st_course_at_this_hospital_desc == "NONE"]
# we have dates and info on treatment given elsewhere

# duplicate ids?
which(duplicated(serous_data$mrn))

serous_data <- ovary_data %>% 
  mutate(mrn = as.character(mrn)) %>% 
  filter(histology_desc %in% c(
    "CARCINOMA NOS", "MALIGNANT SEROUS TUMOR",
    "PAPILLARY SEROUS CYSTADENOCARCINOMA",
    "SEROUS CYSTADENOCARCINOMA NOS", "SEROUS SURFACE PAPILLARY CARCINOMA"
  )) %>% 
  select(mrn, tumor_id, tumor_seq_num,
         cdsc_primary_site_group, histology_desc,
         dx_dt,
         first_treatment_dt, 
         summary_of_rx_1st_course_desc, 
         systemic_rx_surgery_seq_num,
         drug_given_as_part_of_therapy_desc,
         rx_summary_treatment_status_desc)

serous_data1 <-
  serous_data %>% 
  mutate(treatment_count = sapply(strsplit(
    summary_of_rx_1st_course_desc, "/"), length), 
    .before = summary_of_rx_1st_course_desc) %>% 
  separate_wider_delim(cols = summary_of_rx_1st_course_desc, delim = "/",
                       names = c(paste0("treatment_type_sequence_", 1:max(.$treatment_count), names_sep = "")), 
                       too_few = "align_start", too_many = "merge", 
                       cols_remove = FALSE) %>% 
  select(-treatment_count) %>% 
  mutate(first_treatment_type_of_interest = case_when(
    str_detect(treatment_type_sequence_1, "SURG")         ~ "Surgery",
    str_detect(treatment_type_sequence_1, "CHEM")         ~ "Chemotherapy"
  ), .before = treatment_type_sequence_1) %>% 
  rename(first_treatment_dt_including_non_surg_non_chemo = first_treatment_dt) %>% 
  mutate(first_treatment_dt_of_interest = case_when(
    !is.na(first_treatment_type_of_interest)              ~ first_treatment_dt_including_non_surg_non_chemo
  ), .after = first_treatment_type_of_interest) %>% 
  
  mutate(exclusion_no_treatment_of_interest = case_when(
    is.na(first_treatment_type_of_interest)               ~ "excluded",
    TRUE                                                  ~ "include"
  ), .before = first_treatment_type_of_interest) %>% 
  
  group_by(mrn) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  arrange(desc(n), exclusion_no_treatment_of_interest, mrn, tumor_seq_num) %>% 
  mutate(exclusion_no_treatment_of_interest = case_when(
    exclusion_no_treatment_of_interest == "excluded"       ~ "excluded",
    TRUE                                                  ~ NA_character_
  )) %>% 
  
  select(mrn : tumor_seq_num,
         dx_dt, exclusion_no_treatment_of_interest,
         first_treatment_type_of_interest :
           summary_of_rx_1st_course_desc, 
         first_treatment_dt_including_non_surg_non_chemo :
           rx_summary_treatment_status_desc,
         systemic_rx_surgery_seq_num,
         cdsc_primary_site_group,
         histology_desc)


# save

write_rds(serous_data1, 
          paste0("data/", "Sigma_HGSC_patients_data_", today(), ".rds"))
write_csv(serous_data1, 
          paste0("data/", "Sigma_HGSC_patients_data_", today(), ".csv"))

write_rds(serous_data1, 
          paste0(path, "/SIGMA",
                 "/processed data/",
                 "Sigma_HGSC_patients_data_", today(), ".rds"))
write_csv(serous_data1, 
          paste0(path, "/SIGMA",
                 "/processed data/",
                 "Sigma_HGSC_patients_data_", today(), ".csv"))



