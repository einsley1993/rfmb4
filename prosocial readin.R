
################################
#### Sep. 15, 2024 by Yi Li ####
################################

getwd()
setwd("/Users/yili/Desktop/Meta Analysis and Replication Forecasting/MB4")

library(readxl)
# MA dataset:
Prosocial_MA_raw <- read_excel("Prosocial MA.xlsx")
Prosocial_MA_raw$study_ID

# REP dataset:
REP_redacted <- read.csv("~/Desktop/Meta Analysis and Replication Forecasting/MB4/REP_redacted.csv")

REP_redacted <- REP_redacted %>% 
  mutate(vi = d_var_calc) # only variance


library(tidyverse)


# Exclude row 123: study_ID = "Van_de_Vondervoort_Hamlin_2018a", because of 
# the special study design.
# they made a comparison between two conditions (prosocial, antisocial) and the 
# number of children who selected the puppet for each condition, which is why 
# there is both n_1 and n_2.
Prosocial_MA_filt <- Prosocial_MA_raw %>% 
  filter(dependent_measure == "target_selection", 
         study_ID != "Van_de_Vondervoort_Hamlin_2018a")

Prosocial_MA_v1 <- Prosocial_MA_filt %>% 
  mutate(haldane_correction = ifelse((x_1 == n_1) | (x_1 == 0), 0.5, 0)) %>% 
  mutate(n_a = x_1 + haldane_correction,
         n_b = n_1 - x_1 + haldane_correction,
         n_c = x_2 + haldane_correction,
         n_d = n_1 - x_2 + haldane_correction,
         n_2 = n_1) %>% 
  mutate(logOR = log( (n_a*n_d) / (n_b*n_c) ),
         SE.logOR = sqrt( (1/n_a) + (1/n_b) + (1/n_c) + (1/n_d) ),
         d_calc = logOR*sqrt(3)/pi,
         d_var_calc = (SE.logOR^2) * 3/(pi^2)) %>% 
  select(study_ID, d_calc, d_var_calc,
         expt_condition, response_mode, exposure_phase, method,
         dependent_measure, participant_design, native_lang, infant_type, 
         n_1, n_2, mean_age_1, mean_age_2, 
         x_1, x_2, SD_1, SD_2, 
         scenarios, stimuli, stimuli_type, choice_object, Hamlin_Lab, 
         intent_val, outcome_val, n_excluded_1, n_excluded_2) %>% 
  mutate(yi = d_calc,
         vi = d_var_calc) 

# Compute p-value for each person
Prosocial_MA <- Prosocial_MA_v1 %>%
  mutate(a = x_1, b = n_1 - x_1, c = x_2, d = n_2 - x_2) %>% 
  mutate(c = as.integer(c),
         d = as.integer(d)) %>% 
  rowwise() %>%
  mutate(p_value = {
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)
    test_result <- fisher.test(contingency_table)
    test_result$p.value
  }) %>%
  ungroup()


# P-value distribution
hist(Prosocial_MA$p_value)
summary(Prosocial_MA$p_value)

count_significant <- Prosocial_MA %>%
  filter(p_value < 0.05) %>%
  summarise(count = n()) %>%
  pull(count)



hist(Prosocial_MA$yi)  
hist(Prosocial_MA$vi) 
summary(Prosocial_MA$vi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02151 0.11028 0.15705 0.17895 0.20069 0.73568 

# compare d_var_calc to REP_redacted:
hist(REP_redacted$vi)
summary(REP_redacted$vi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1105  0.2584  0.3546  0.4975  0.6206  1.6211 

write.csv(Prosocial_MA, "mb4_dma.csv")
write.csv(REP_redacted, "mb4_dr.csv")




