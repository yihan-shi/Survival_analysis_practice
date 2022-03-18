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




# censor time set to 3 years ----------------------------------------------
# break down proportion of revision patients
# check gbp & dds
# try 1) get rid of dds
# 2) artifically censor everyone at 3 years after index
# 3) get rid of anyone who was censored before three years (but if they had
#    an event before 3 years, but didn't have enough follow up after that,
#    thats ok too) - so if they had event at 1 year but only had 2 years follow up, keep them
# 3) restrict patients to follow up < 3 years (en_rlast - index <= 3 years)

# do the same for 2 years


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
         status = 0) %>%
  # these patients only have 1 outcome, so artificially censor at 3 years
  mutate(censor_time = ifelse(time < 1096, NA, 1096)) %>%
  filter(!is.na(censor_time))


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

## slg = index surgery
slg <- all_patients %>%
  filter(occurrence == 1 & type == "slg")
slg_patient <- slg$PAT_ID

## agb = index surgery
agb <- all_patients %>%
  filter(occurrence == 1 & type == "agb")
agb_patient <- agb$PAT_ID

## dds = index surgery
dds <- all_patients %>%
  filter(occurrence == 1 & type == "dds")
dds_patient <- dds$PAT_ID
dds_surv <- all_patients %>%
  filter(PAT_ID %in% dds_patient)


## combine & limit to the first event after index surgery
# km_all <- rbind(dds_surv, agb_surv, slg_surv, gbp_surv)
km_all <- rbind(dds, agb, slg, gbp)
view(km_all)

# Kaplan-Meier -------------------------------------------------------
# Produce the Kaplan-Meier estimates of the probability of survival
# time to the first event

km <- with(km_all, Surv(time, status) ~ type)
km_fit <- survfit(Surv(time, status) ~ type, data = km_all)

# Print the estimates for 1, 30, 60 and 90 days, and then every 90 days thereafter
summary(km_fit, times = c(1, 300, 800, 1000*(1:4)))

par(mar = c(3,3,3,3))
plot(km_fit,
     main = 'Kaplan Meyer Plot for bariatric surgery revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"),
     lwd = 2)
legend("topright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))



# Nelson-Aalen ------------------------------------------------------------
na_fit <- survfit(Surv(time, status) ~ type, type = "fh", data = km_all)
summary(na_fit, times = c(1, 300, 800, 1000*(1:4)))

par(mar = c(3,3,3,3))
plot(na_fit,
     main = 'Nelson-Aalen Plot for bariatric surgery revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"),
     lwd = 2)
legend("topright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))

# Cox exp---------------------------------------------------------------------
exp_plot <- survreg(Surv(time, status) ~ type,
                    data = km_all, dist = "exponential")
plot(predict(exp_plot,
             newdata = list(type = "agb"),
             type = "quantile",
             p = seq(.01,.99, by=.01)),
     seq(.99,.01,by = -.01),
     col = "#000000",
     type = "l",
     main = "Cox proportional model based on exponential distribution",
     xlab = "Time",
     ylab = "Estimated survivor function",
     lty = "solid",
     lwd = 2)
lines(predict(exp_plot,
              newdata = list(type = "dds"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#E69F00",
      lty = "dashed")
lines(predict(exp_plot,
              newdata = list(type = "gbp"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#56B4E9",
      lty = "dotted")
lines(predict(exp_plot,
              newdata = list(type = "slg"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#009E73",
      lty = "dotdash")
legend("topright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))


# Cox weibull -----------------------------------------------------------------
weib_plot <- survreg(Surv(time, status) ~ type,
                     data = km_all, dist = "weibull")
plot(predict(weib_plot,
             newdata = list(type = "agb"),
             type = "quantile",
             p = seq(.01,.99, by=.01)),
     seq(.99,.01,by = -.01),
     col = "#000000",
     type = "l",
     main = "Cox proportional model based on exponential distribution",
     xlab = "Time",
     ylab = "Estimated survivor function",
     lty = "solid",
     lwd = 2)
lines(predict(weib_plot,
              newdata = list(type = "dds"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#E69F00",
      lty = "dashed")
lines(predict(weib_plot,
              newdata = list(type = "gbp"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#56B4E9",
      lty = "dotted")
lines(predict(weib_plot,
              newdata = list(type = "slg"),
              type="quantile",
              p = seq(.01,.99,by = .01)),
      seq(.99,.01,by  =-.01),
      col = "#009E73",
      lty = "dotdash")
legend("topright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))


# cox fit -----------------------------------------------------------------
ggsurvplot(km_fit, data = km_all, risk.table = TRUE, size = 1,
           palette = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
           legend = "top",
           linetype = "strata",
           fun = "pct",
           legend.title = "Surgery type",
           legend.labs = c("agb", "dds", "gbp", "slg"),
           risk.table.height = 0.32)

# cum hazard
ggsurvplot(km_fit, data = km_all, fun = "cumhaz",
           palette = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
           legend.title = "Surgery type")

# table summary
coxph(Surv(time, status) ~ type, data = km_all) %>%
  gtsummary::tbl_regression(exp = TRUE)


# recurrent event model ---------------------------------------------------------
# patients who experienced recurrent events

# patient id who have revisions after first index surgeries: multiple_patient
recurrent <- surgeries %>%
  filter(PAT_ID %in% multiple_patient) %>%
  group_by(PAT_ID) %>%
  pivot_longer(
    cols = c("gbp", "slg", "agb", "dds", "gjs",
             "rev_of_agb", "rem_of_agb", "hhr", "ejs", "rev_of_gjs"),
    names_to = "type") %>%
  filter(value == 1) %>%
  select(-value) %>%
  mutate(time = as.numeric(FROM_DT - lag(index)),
         status = 1) %>% # observed event
  filter(!time %in% c(NA, 0))

# combine patient types
rec_all <- rbind(uniq_cases, recurrent) %>%
  group_by(PAT_ID) %>%
  mutate(occurrence = row_number()) %>%
  filter(time > 0)

view(rec_all)


## patients who experienced recurrent events
p_recurrent <- rec_all %>%
  group_by(PAT_ID) %>%
  filter(status == 1 & n() > 1) %>%
  mutate(diff = time - lag(time)) %>%
  mutate(diff = ifelse(status == 1 & is.na(diff), 0, diff))

# create a list of unique patient ids
patients <- unique(p_recurrent$PAT_ID)

# create empty dataframe
output_rec <- data.frame(matrix(ncol = 8, nrow = 0))
x <- c("PAT_ID", "time", "diff", "n", "tstart", "tstop", "enr_last", "index")
colnames(output_rec) <- x


for(i in 1:length(unique(patients))){
  # create a separate df for 1 patient
  df <- p_recurrent %>%
    filter(PAT_ID == unique(patients)[i])

  df <- df %>%
    mutate(n = 1:n()) %>%
    mutate(tstop = ifelse(n == 1, df[n+1,]$diff - 1, 0),
           tstart = 0) %>%
    select("PAT_ID", "time", "diff", "n", "tstart", "tstop", "enr_last", "index")

  for(j in 2:dim(df)[1]){
    df$tstart[j] <- df[j-1,]$tstop + 1
    df$tstop[j] <- df[j,]$tstart + df[j+1,]$diff - 1
  }

  df$tstop[dim(df)[1]] <- as.numeric(difftime(df$enr_last[dim(df)[1]],
                                              df$index[dim(df)[1]],
                                              units = c("days")))
  output_rec <- rbind(output_rec, df)
}

output_rec <- output_rec %>%
  mutate(tstop = ifelse(tstop < 0, 0, tstop)) %>%
  select("PAT_ID", "n", "tstart", "tstop")
view(output_rec)


## patients who experienced 1 event
p_nonrecurrent <-  rec_all %>%
  group_by(PAT_ID) %>%
  filter(status == 1 & n() == 1) %>%
  mutate(tstart = 0,
         tstop = as.numeric(difftime(FROM_DT, index, units = c("days"))),
         n = 1) %>%
  select("PAT_ID", "n", "tstart", "tstop")

## censored patient (no event)
p_no <- rec_all %>%
  group_by(PAT_ID) %>%
  filter(status == 0) %>%
  mutate(tstart = 0,
         tstop = as.numeric(difftime(enr_last, index, units = c("days"))),
         n = 1) %>%
  select("PAT_ID", "n", "tstart", "tstop")

output_all <- rbind(output_rec, p_nonrecurrent, p_no) %>%
  filter(tstop > 0 & (tstart < tstop))


# add surgery type --------------------------------------------------------
recur_data <- output_all %>%
  inner_join(rec_all, c("PAT_ID" = "PAT_ID")) %>%
  distinct(PAT_ID, tstart, tstop, status, type, .keep_all = TRUE) %>%
  group_by(PAT_ID)


# Andersen-Gill (AG) model ------------------------------------------------
model_AG <- coxph(Surv(tstart, tstop, status) ~ type + cluster(PAT_ID),
                  method = "breslow",
                  cluster = PAT_ID,
                  robust = TRUE, data = recur_data)
summary(model_AG)

# plot AG
plot(survfit(model_AG),
     main = 'Andersen-Gill Plot for recurrent revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"),
     lwd = 1)
legend("bottomright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))

# table summary of AG
model_AG %>%
  gtsummary::tbl_regression(exp = TRUE)


# PWP-TT model ----------------------------------------------------------------------
model_pwptt <- coxph(Surv(tstart, tstop, status) ~ type + cluster(PAT_ID) + strata(n),
                     method = "breslow",
                     data = recur_data)
summary(model_pwptt)

# plot PWP-TT
plot(survfit(model_pwptt),
     main = 'PWP-TT Plot for recurrent revisions',
     xlab = "Revision time",
     ylab = "Survival Probabilities",
     col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
     lty = c("solid", "dashed", "dotted", "dotdash"),
     lwd = 1)
legend("bottomright",
       c("agb", "dds", "gbp", "slg"),
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73"),
       lty = c("solid", "dashed", "dotted", "dotdash"))


model_pwptt %>%
  gtsummary::tbl_regression(exp = TRUE)

