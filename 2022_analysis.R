# look for differences between 4 types of patients
# separated by which type of surgery was their “index surgery”:

  # Gbp (gastric bypass)
  # Slg (sleeve gastrectomy)
  # Agb (adjustable gastric band)
  # Dds (duodenal switch)

# see whether they have these any type of events afterwards on a different date
# from their index surgery

# Gjs (gastrojejunostomy)
# Rev_of_agb (revision of agb)
# Rem_of_agb (removal of agb – do not count this if it is the ONLY code after the initial agb and not accompanied by something else, either the same day or later)
# Hhr (hiatal hernia repair)
# Ejs (esophagojejunostomy)
# Rev_of_gjs (revision of gjs)

# load libraries ----------------------------------------------------------

library(survival)
library(survminer)
library(ranger)
library(gt)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggfortify)
library(devtools)
# library(OIsurv)
library(lubridate)


# load data ---------------------------------------------------------------

# data <- load("Desktop/Duke_2021-22/Jiang research/Survival_analysis_practice/data/revisions.RData")
data <- load("data/revisions.RData")
view(repeats)
view(enroll)
view(events)

dim(repeats) # 7522 rows
dim(enroll) # 48997 rows
dim(events) # 53143 rows

# clean data --------------------------------------------------------------
surgeries <- events %>%
  inner_join(enroll, c("PAT_ID" = "pat_id")) %>%
  select("PAT_ID", "index", "FROM_DT", "gbp", "slg", "agb", "dds", "gjs",
         "rev_of_agb", "rem_of_agb", "hhr", "ejs", "rev_of_gjs", "enr_last")


# revision patients (non-recurrent) -------------------------------------------------------
multiple <- surgeries %>%
  group_by(PAT_ID) %>%
  count() %>%
  filter(n > 1)

# patient who have revisions after first index surgeries
multiple_patient <- multiple$PAT_ID

repeated <- surgeries %>%
  filter(PAT_ID %in% multiple_patient) %>%
  group_by(PAT_ID) %>%
  slice(1:3) %>% # choose the first 2 events for non-recurrent event analysis
  pivot_longer(
    cols = c("gbp", "slg", "agb", "dds", "gjs",
             "rev_of_agb", "rem_of_agb", "hhr", "ejs", "rev_of_gjs"),
    names_to = "type") %>%
  filter(value == 1) %>%
  select(-value) %>%
  mutate(time = as.numeric(FROM_DT - lag(index)),
         status = 1) %>% # observed event
  filter(!time %in% c(NA, 0))

# view(repeated)

# no revision patients ----------------------------------------------------
uniq <- surgeries %>%
  group_by(PAT_ID) %>%
  count() %>%
  filter(n == 1)

# patient who doesn't have revisions after first index surgeries
uniq_patient <- uniq$PAT_ID

uniq_cases <- surgeries %>%
  filter(PAT_ID %in% uniq_patient) %>%
  group_by(PAT_ID) %>%
  pivot_longer(
    cols = c("gbp", "slg", "agb", "dds", "gjs",
             "rev_of_agb", "rem_of_agb", "hhr", "ejs", "rev_of_gjs"),
    names_to = "type") %>%
  filter(value == 1) %>%
  select(-value) %>%
  mutate(time = as.numeric(difftime(enr_last, index, units = c("days"))),
         status = 0) # censored


# combine patient types ---------------------------------------------------
all_patients <- rbind(uniq_cases, repeated) %>%
  group_by(PAT_ID) %>%
  mutate(occurrence = row_number()) %>%
  filter(time > 0)

view(all_patients)

## gbp = index surgery
gbp <- all_patients %>%
  filter(occurrence == 1 & type == "gbp")
gbp_patient <- gbp$PAT_ID
gbp_surv <- all_patients %>%
  filter(PAT_ID %in% gbp_patient)

## slg = index surgery
slg <- all_patients %>%
  filter(occurrence == 1 & type == "slg")
slg_patient <- slg$PAT_ID
slg_surv <- all_patients %>%
  filter(PAT_ID %in% slg_patient)

## agb = index surgery
agb <- all_patients %>%
  filter(occurrence == 1 & type == "agb")
agb_patient <- agb$PAT_ID
agb_surv <- all_patients %>%
  filter(PAT_ID %in% agb_patient)

## dds = index surgery
dds <- all_patients %>%
  filter(occurrence == 1 & type == "dds")
dds_patient <- dds$PAT_ID
dds_surv <- all_patients %>%
  filter(PAT_ID %in% dds_patient)

## combine

km_all <- rbind(dds, agb, slg, gbp)


# Kaplan-Meier -------------------------------------------------------
# Produce the Kaplan-Meier estimates of the probability of survival
# time to the first event

km <- with(km_all, Surv(time, status) ~ type)
km_fit <- survfit(Surv(time, status) ~ type, data = km_all)

# Print the estimates for 1, 30, 60 and 90 days, and then every 90 days thereafter
summary(km_fit, times = c(1, 300, 800, 1000*(1:4)))

plot(km_fit,
     main = 'Kaplan Meyer Plot for bariatric surgery revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"))
legend("bottomright",
       c("agb", "gbp", "slg", "dds"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))



# Nelson-Aalen ------------------------------------------------------------
na_fit <- survfit(Surv(time, status) ~ type, type = "fh", data = km_all)
summary(na_fit, times = c(1, 300, 800, 1000*(1:4)))

# par(mar=c(1,1,1,1))
plot(na_fit,
     main = 'Nelson-Aalen Plot for bariatric surgery revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"))
legend("bottomright",
       c("agb", "gbp", "slg", "dds"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))

# Cox ---------------------------------------------------------------------







