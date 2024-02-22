library(data.table)
library(tidyverse)
library(dplyr)
library(survial)
library(epiR)

setwd("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB")
dat.baseline <- fread("ukbb_RA&CAD_baselibe.csv")
dat.followup <- fread("ukbb_RA&CAD_followup_2.csv")

dat <- merge(dat.baseline, dat.followup, by.x = "ID", by.y = "f.eid")#502407
rm(dat.baseline)
rm(dat.followup)

# outcome: CVD exposure:RA
dat$followup_event <- ifelse(dat$Stroke == 1 | dat$CAD == 1 | dat$AF == 1 | dat$HF == 1, 1, 0)
table(dat$followup_event)
#    0      1 
#415493  86914 

dat$CVD_difftime <- apply(dat[,c("Stroke_difftime","CAD_difftime","AF_difftime","HF_difftime")], 1, min, na.rm = T)
dat$followup_time <- apply(dat[, c("losetime_difftime","endtime_difftime","CVD_difftime")],1,min, na.rm=T)                                      



# Step 1: only include European individuals
dat <- subset(dat, Enthnic=="1001"|Enthnic=="1002"|Enthnic=="1003"|Enthnic=="1")#472609

# Step 2: exclude individuals with prevalent CVD (stroke, HF, AF, CAD) at baseline
dat <- subset(dat, dat$CVD_difftime > 0)#442657
dat$CVD_selfreport <- ifelse(dat$Stroke_selfreport == 1 | dat$CAD_selfreport == 1 | dat$AF_selfreport == 1 | dat$HF_selfreport == 1, 1,0)
dat <- subset(dat, dat$CVD_selfreport == 0) #441096 some self-reported cases did not provide time, exclude them


dat$exposure_RA <- ifelse(dat$RA ==1 & dat$RA_difftime < dat$CVD_difftime , 1, 0)
table(dat$exposure_RA)

#################################################################################
############# incidence density #######
table(dat$exposure_RA)
table(dat$exposure_RA, by = dat$followup_event)

dat.ra <- subset(dat, dat$exposure_RA == 1)
dat.nonra <- subset(dat, dat$exposure_RA == 0)


#RA
ncase <- 1467
ntar <- sum(dat.ra$followup_time)/365
tmp <- as.matrix(cbind(ncase, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 9041, design = 1,  conf.level = 0.95)*1000

#nonRA
ncase <- 47851
ntar <- sum(dat.nonra$followup_time)/365
tmp <- as.matrix(cbind(ncase, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 432055, design = 1,  conf.level = 0.95)*1000

#################################################################################
############### Basic characteristics #######
library(pastecs)
library(doBy)
library(psych)
library(simputation)
library(mice)
library(janitor)
library(Hmisc)
library(dbplyr)
# basic characteristics
# age
describeBy(dat$Age,list(dat$exposure_RA))
t.test(dat$Age ~ dat$exposure_RA)

#sex
table(dat$sex, by = dat$exposure_RA)
re.dat <- dat %>%
  tabyl(exposure_RA, sex) %>%
  adorn_percentages("row") %>% 
  adorn_pct_formatting(digits = 1)  %>%
  adorn_ns()
re.dat
chisq.test(dat$sex, dat$exposure_RA)

#BMI
dat$BMI_group <- ifelse(dat$BMI >30, 4,
                        ifelse(dat$BMI >25 & dat$BMI <=29.9,3,
                               ifelse(dat$BMI > 18.5 & dat$BMI <=24.9,2,1)))

dat$BMI_group[is.na(dat$BMI_group)] <- 9
re.dat <- dat %>%
  tabyl(exposure_RA, BMI_group) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$BMI_group, dat$exposure_RA)

# smoking
dat$Smoking[dat$Smoking== -3] <- 9
re.dat <- dat %>%
  tabyl(exposure_RA, Smoking) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Smoking, dat$exposure_RA)

#alcohol consumption#Ref1:Never
dat <- dat %>% mutate(Alcoholgroup=if_else(Alcohol=="6",1,
                                           ifelse(Alcohol=="5",2,
                                                  ifelse(Alcohol=="4",3,
                                                         ifelse(Alcohol=="3",4,
                                                                ifelse(Alcohol=="2",5,
                                                                       ifelse(Alcohol=="1",6, 9)))))))


re.dat <- dat %>%
  tabyl(exposure_RA, Alcoholgroup) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Alcoholgroup, dat$exposure_RA)



#Education
dat <- dat %>% mutate(educationgroup=ifelse(Education.qua == "1" | Education.qua == "6",3,
                                                ifelse(Education.qua == "2"| Education.qua == "3"| Education.qua == "4"|Education.qua == "5",2,
                                                       ifelse(Education.qua == "-7",1, 9))))


dat$educationgroup[is.na(dat$educationgroup)] <- 9
re.dat <- dat %>%
  tabyl(exposure_RA, educationgroup) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$educationgroup, dat$exposure_RA)




# hypertension
dat$Hyp_baseline <- ifelse(dat$Hyp == 1 & dat$Hyp_difftime < 0 , 1,0) # baseline Hypertension

re.dat <- dat %>%
  tabyl(exposure_RA, Hyp_baseline) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Hyp_baseline, dat$exposure_RA)




#T2D
dat$T2D_baseline <- ifelse(dat$T2D == 1 & dat$T2D_difftime < 0 , 1,0) # baseline Hypertension

re.dat <- dat %>%
  tabyl(exposure_RA, T2D_baseline) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$T2D_baseline, dat$exposure_RA)

# METs
describeBy(dat$METs,list(dat$exposure_RA))

dat %>%
  group_by(exposure_RA) %>%
  summarise(quant = quantile(METs, probs = c(0.25,0.75)))
wilcox.test(dat$METs, dat$exposure_RA)

# TC
dat$logTC <- log(dat$TC, 2)
describeBy(dat$logTC,list(dat$exposure_RA))
t.test(dat$logTC ~ dat$exposure_RA)


################################################################################
############## Association analysis ##################################
library(survival)
### RA and incident CVD
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF
dat$followup_time_AF <- apply(dat[, c("losetime_difftime","endtime_difftime","AF_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat$followup_time_AF)
status <- as.numeric(dat$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])
### RA and incident CAD
dat$followup_time_CAD <- apply(dat[, c("losetime_difftime","endtime_difftime","CAD_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat$followup_time_CAD)
status <- as.numeric(dat$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident HF
dat$followup_time_HF <- apply(dat[, c("losetime_difftime","endtime_difftime","HF_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat$followup_time_HF)
status <- as.numeric(dat$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
dat$followup_time_Stroke <- apply(dat[, c("losetime_difftime","endtime_difftime","Stroke_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat$followup_time_Stroke)
status <- as.numeric(dat$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_RA + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2

re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])
re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)

################################################################################
#################### Sensitivity analyses ##############################
# sex specific
dat.sex <- subset(dat, dat$sex == 0) # female 247794
dat.sex <- subset(dat, dat$sex == 1) # male
### RA and incident CVD
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$Hyp_baseline + dat.sex$T2D_baseline + dat.sex$METs + dat.sex$TC, data=dat.sex))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF
                                      
time <-  as.numeric(dat.sex$followup_time_AF)
status <- as.numeric(dat.sex$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$Hyp_baseline + dat.sex$T2D_baseline + dat.sex$METs + dat.sex$TC, data=dat.sex))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident CAD
                                      
time <-  as.numeric(dat.sex$followup_time_CAD)
status <- as.numeric(dat.sex$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$Hyp_baseline + dat.sex$T2D_baseline + dat.sex$METs + dat.sex$TC, data=dat.sex))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident HF
                                     
time <-  as.numeric(dat.sex$followup_time_HF)
status <- as.numeric(dat.sex$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$Hyp_baseline + dat.sex$T2D_baseline + dat.sex$METs + dat.sex$TC, data=dat.sex))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident stroke
                                    
time <-  as.numeric(dat.sex$followup_time_Stroke)
status <- as.numeric(dat.sex$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_RA + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$Hyp_baseline + dat.sex$T2D_baseline + dat.sex$METs + dat.sex$TC, data=dat.sex))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)
rm(dat.sex)
# age group based on 60
dat.age <- subset(dat, dat$Age < 60) #253579
dat.age <- subset(dat, dat$Age >= 60)

### RA and incident CVD
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$Hyp_baseline + dat.age$T2D_baseline + dat.age$METs + dat.age$TC, data=dat.age))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF
dat.age$followup_time_AF <- apply(dat.age[, c("losetime_difftime","endtime_difftime","AF_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat.age$followup_time_AF)
status <- as.numeric(dat.age$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$Hyp_baseline + dat.age$T2D_baseline + dat.age$METs + dat.age$TC, data=dat.age))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident CAD
dat.age$followup_time_CAD <- apply(dat.age[, c("losetime_difftime","endtime_difftime","CAD_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat.age$followup_time_CAD)
status <- as.numeric(dat.age$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$Hyp_baseline + dat.age$T2D_baseline + dat.age$METs + dat.age$TC, data=dat.age))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident HF
dat.age$followup_time_HF <- apply(dat.age[, c("losetime_difftime","endtime_difftime","HF_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat.age$followup_time_HF)
status <- as.numeric(dat.age$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$Hyp_baseline + dat.age$T2D_baseline + dat.age$METs + dat.age$TC, data=dat.age))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident stroke
dat.age$followup_time_Stroke <- apply(dat.age[, c("losetime_difftime","endtime_difftime","Stroke_difftime")],1,min, na.rm=T)                                      
time <-  as.numeric(dat.age$followup_time_Stroke)
status <- as.numeric(dat.age$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_RA + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$Hyp_baseline + dat.age$T2D_baseline + dat.age$METs + dat.age$TC, data=dat.age))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)
rm(dat.age)

# exclude the first one years
dat.ex1 <- dat[!(which(dat$followup_event == 1 & dat$followup_time <365)),] #438550

### RA and incident CVD
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$Hyp_baseline + dat.ex1$T2D_baseline + dat.ex1$METs + dat.ex1$TC, data=dat.ex1))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident AF
dat.ex1 <- dat[!(which(dat$AF == 1 & dat$followup_time_AF <365)),] #

time <-  as.numeric(dat.ex1$followup_time_AF)
status <- as.numeric(dat.ex1$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$Hyp_baseline + dat.ex1$T2D_baseline + dat.ex1$METs + dat.ex1$TC, data=dat.ex1))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident CAD
dat.ex1 <- dat[!(which(dat$CAD == 1 & dat$followup_time_CAD <365)),] #4

time <-  as.numeric(dat.ex1$followup_time_CAD)
status <- as.numeric(dat.ex1$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$Hyp_baseline + dat.ex1$T2D_baseline + dat.ex1$METs + dat.ex1$TC, data=dat.ex1))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident HF
dat.ex1 <- dat[!(which(dat$HF == 1 & dat$followup_time_HF <365)),] #4

time <-  as.numeric(dat.ex1$followup_time_HF)
status <- as.numeric(dat.ex1$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$Hyp_baseline + dat.ex1$T2D_baseline + dat.ex1$METs + dat.ex1$TC, data=dat.ex1))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
dat.ex1 <- dat[!(which(dat$Stroke == 1 & dat$followup_time_Stroke <365)),] #4

time <-  as.numeric(dat.ex1$followup_time_Stroke)
status <- as.numeric(dat.ex1$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_RA + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$Hyp_baseline + dat.ex1$T2D_baseline + dat.ex1$METs + dat.ex1$TC, data=dat.ex1))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)

rm(dat.ex1)

#exclude kindship
dat.exkin <- subset(dat, dat$kinship == 0) #295159

### RA and incident CVD
time <-  as.numeric(dat.exkin$followup_time)
status <- as.numeric(dat.exkin$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$Hyp_baseline + dat.exkin$T2D_baseline + dat.exkin$METs + dat.exkin$TC, data=dat.exkin))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident AF
dat.exkin<- dat[!(which(dat$AF == 1 & dat$followup_time_AF <365)),] #

time <-  as.numeric(dat.exkin$followup_time_AF)
status <- as.numeric(dat.exkin$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$Hyp_baseline + dat.exkin$T2D_baseline + dat.exkin$METs + dat.exkin$TC, data=dat.exkin))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident CAD
dat.exkin <- dat[!(which(dat$CAD == 1 & dat$followup_time_CAD <365)),] #4

time <-  as.numeric(dat.exkin$followup_time_CAD)
status <- as.numeric(dat.exkin$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$Hyp_baseline + dat.exkin$T2D_baseline + dat.exkin$METs + dat.exkin$TC, data=dat.exkin))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident HF
dat.exkin <- dat[!(which(dat$HF == 1 & dat$followup_time_HF <365)),] #4

time <-  as.numeric(dat.exkin$followup_time_HF)
status <- as.numeric(dat.exkin$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$Hyp_baseline + dat.exkin$T2D_baseline + dat.exkin$METs + dat.exkin$TC, data=dat.exkin))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
dat.exkin <- dat[!(which(dat$Stroke == 1 & dat$followup_time_Stroke <365)),] #4

time <-  as.numeric(dat.exkin$followup_time_Stroke)
status <- as.numeric(dat.exkin$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_RA + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$Hyp_baseline + dat.exkin$T2D_baseline + dat.exkin$METs + dat.exkin$TC, data=dat.exkin))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)



# propensity score
library(MatchIt) 
library(survival)
dat.propensity <- dat[,c("ID","sex","Age","Townsend depreivation index","exposure_RA","followup_event")]#441096
dat.propensity.nomissing <- dat.propensity[complete.cases(dat.propensity), ] #440569

dat.propensity.match <- matchit(exposure_RA ~ sex + Age + `Townsend depreivation index`, data = dat.propensity.nomissing,  distance = "glm", method = "nearest", ratio = 4)
dat.propensity.match.data <- match.data(dat.propensity.match)#45135

dat.propensity.match.data <- as.data.frame(dat.propensity.match.data$ID)

dat.propensity.match.data <- merge(dat.propensity.match.data, dat, by.x = "dat.propensity.match.data$ID", by.y = "ID")

### RA and incident CVD
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$Hyp_baseline + dat.propensity.match.data$T2D_baseline + dat.propensity.match.data$METs + dat.propensity.match.data$TC, data=dat.propensity.match.data))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF

time <-  as.numeric(dat.propensity.match.data$followup_time_AF)
status <- as.numeric(dat.propensity.match.data$AF)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$Hyp_baseline + dat.propensity.match.data$T2D_baseline + dat.propensity.match.data$METs + dat.propensity.match.data$TC, data=dat.propensity.match.data))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident CAD

time <-  as.numeric(dat.propensity.match.data$followup_time_CAD)
status <- as.numeric(dat.propensity.match.data$CAD)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$Hyp_baseline + dat.propensity.match.data$T2D_baseline + dat.propensity.match.data$METs + dat.propensity.match.data$TC, data=dat.propensity.match.data))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident HF

time <-  as.numeric(dat.propensity.match.data$followup_time_HF)
status <- as.numeric(dat.propensity.match.data$HF)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$Hyp_baseline + dat.propensity.match.data$T2D_baseline + dat.propensity.match.data$METs + dat.propensity.match.data$TC, data=dat.propensity.match.data))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident stroke

time <-  as.numeric(dat.propensity.match.data$followup_time_Stroke)
status <- as.numeric(dat.propensity.match.data$Stroke)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_RA + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$Hyp_baseline + dat.propensity.match.data$T2D_baseline + dat.propensity.match.data$METs + dat.propensity.match.data$TC, data=dat.propensity.match.data))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)


