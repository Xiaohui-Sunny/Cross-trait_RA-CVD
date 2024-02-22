library(data.table)
library(tidyverse)
library(dplyr)
library(survial)

setwd("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB")
dat.baseline <- fread("ukbb_RA&CAD_baselibe.csv")
dat.followup <- fread("ukbb_RA&CAD_followup_2.csv")

dat <- merge(dat.baseline, dat.followup, by.x = "ID", by.y = "f.eid")#502407
rm(dat.baseline)
rm(dat.followup)


# outcome: RA exposure:CVD

dat$followup_event <- ifelse(dat$RA == 1, 1, 0)
table(dat$followup_event)
#  0      1 
#490577  11830

dat$CVD_difftime <- apply(dat[,c("Stroke_difftime","CAD_difftime","AF_difftime","HF_difftime")], 1, min, na.rm = T)
dat$CVD <- ifelse(dat$Stroke == 1 | dat$CAD == 1 | dat$HF == 1 | dat$AF == 1 ,1, 0)
dat$followup_time <- apply(dat[, c("losetime_difftime","endtime_difftime","RA_difftime")],1,min, na.rm=T)                                      



# Step 1: only include European individuals
dat <- subset(dat, Enthnic=="1001"|Enthnic=="1002"|Enthnic=="1003"|Enthnic=="1")#472609

# Step 2: exclude individuals with prevalent RA at baseline
dat <- subset(dat, dat$RA_difftime > 0)#466527
dat <- subset(dat, dat$RA_selfreport == 0) #466285 some self-reported cases did not provide time, exclude them


# Step 3: exclude the individuals first diagnosis of RA and subsequently diagnosis of CVD
#dat <- subset(dat, !(dat$CVD ==1 & dat$RA == 1 & dat$RA_difftime < dat$CVD_difftime)) #465821
#dat$exposure_CVD <- ifelse(dat$CVD == 1,1, 0)

dat$exposure_CVD <- ifelse(dat$CVD == 1 & dat$CVD_difftime < dat$RA_difftime , 1, 0)

#################################################################################
############# incidence density #######
table(dat$exposure_CVD)
table(dat$exposure_CVD, by= dat$RA)

dat.CVD <- subset(dat, dat$exposure_CVD == 1) #78282
dat.nonCVD <- subset(dat, dat$exposure_CVD == 0) #388003


#CVD
ncase <- 1012
ntar <- sum(dat.CVD$followup_time)/365
tmp <- as.matrix(cbind(ncase, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 78282, design = 1,  conf.level = 0.95)*1000

#nonCVD
ncase <- 3707
ntar <- sum(dat.nonCVD$followup_time)/365
tmp <- as.matrix(cbind(ncase, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 388003, design = 1,  conf.level = 0.95)*1000


#################################################################################
############### Basic characteristics #######
library(pastecs)
library(doBy)
library(psych)
library(janitor)
library(Hmisc)
library(dbplyr)
# basic characteristics
# age
describeBy(dat$Age,list(dat$exposure_CVD))
t.test(dat$Age ~ dat$exposure_CVD)

#sex
table(dat$sex, by = dat$exposure_CVD)
re.dat <- dat %>%
  tabyl(exposure_CVD, sex) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$sex, dat$exposure_CVD)

#BMI
describeBy(dat$BMI,list(dat$exposure_CVD))

dat$BMI_group <- ifelse(dat$BMI >30, 4,
                        ifelse(dat$BMI >25 & dat$BMI <=29.9,3,
                              ifelse(dat$BMI > 18.5 & dat$BMI <=24.9,2,1)))

dat$BMI_group[is.na(dat$BMI_group)] <- 9
re.dat <- dat %>%
  tabyl(exposure_CVD, BMI_group) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$BMI_group, dat$exposure_CVD)


# smoking
dat$Smoking[dat$Smoking == -3] <- 9
re.dat <- dat %>%
  tabyl(exposure_CVD, Smoking) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Smoking, dat$exposure_CVD)

#alcohol consumption#Ref1:Never
dat <- dat %>% mutate(Alcoholgroup=if_else(Alcohol=="6",1,
                                           ifelse(Alcohol=="5",2,
                                                  ifelse(Alcohol=="4",3,
                                                         ifelse(Alcohol=="3",4,
                                                                ifelse(Alcohol=="2",5,
                                                                       ifelse(Alcohol=="1",6, 9)))))))


re.dat <- dat %>%
  tabyl(exposure_CVD, Alcoholgroup) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Alcoholgroup, dat$exposure_CVD)



#Education
dat <- dat %>% mutate(educationgroup=ifelse(Education.qua == "1" | Education.qua == "6",3,
                                            ifelse(Education.qua == "2"| Education.qua == "3"| Education.qua == "4"|Education.qua == "5",2,
                                                   ifelse(Education.qua == "-7",1, 9))))

dat$educationgroup[is.na(dat$educationgroup)] <- 9
re.dat <- dat %>%
  tabyl(exposure_CVD, educationgroup) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$educationgroup, dat$exposure_CVD)

# hypertension
dat$Hyp_baseline <- ifelse(dat$Hyp == 1 & dat$Hyp_difftime < 0 , 1,0) # baseline Hypertension

re.dat <- dat %>%
  tabyl(exposure_CVD, Hyp_baseline) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$Hyp_baseline, dat$exposure_CVD)

#T2D
dat$T2D_baseline <- ifelse(dat$T2D == 1 & dat$T2D_difftime < 0 , 1,0) # baseline Hypertension

re.dat <- dat %>%
  tabyl(exposure_CVD, T2D_baseline) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 1) %>%
  adorn_ns()
re.dat
chisq.test(dat$T2D_baseline, dat$exposure_CVD)

# METs
describeBy(dat$METs,list(dat$exposure_CVD))

dat %>%
  group_by(exposure_CVD) %>%
  summarise(quant = quantile(METs, probs = c(0.25,0.75)))
wilcox.test(dat$METs, dat$exposure_CVD)

# TC
dat$logTC <- log(dat$TC, 2)
describeBy(dat$logTC,list(dat$exposure_CVD))
t.test(dat$logTC ~ dat$exposure_CVD)

###############################################################################
############## Association analysis ##################################
library(survival)
### CVD and incident RA
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_CVD + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_CVD + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$METs, data=dat))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

###  AF and incident RA
dat$exposure_AF <- ifelse(dat$AF == 1 & dat$AF_difftime < dat$RA_difftime , 1, 0)                                  
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_AF + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_AF + dat$Age + dat$sex + dat$BMI_group  + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$METs, data=dat))
re2

re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident CAD
dat$exposure_CAD <- ifelse(dat$CAD == 1 & dat$CAD_difftime < dat$RA_difftime , 1, 0)                                  
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_CAD + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_CAD + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$METs, data=dat))
re2

re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident HF
dat$exposure_HF <- ifelse(dat$HF == 1 & dat$HF_difftime < dat$RA_difftime , 1, 0)                                  
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_HF + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_HF + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$METs, data=dat))
re2

re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
dat$exposure_Stroke <- ifelse(dat$Stroke == 1 & dat$Stroke_difftime < dat$RA_difftime , 1, 0)                                  
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat$exposure_Stroke + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$exposure_Stroke + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$METs, data=dat))
re2

re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])
re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)

################################################################################
#################### Sensitivity analyses ##############################
# sex specific
dat.sex <- subset(dat, dat$sex == 0) # female 253074
dat.sex <- subset(dat, dat$sex == 1) # male
### CVD and incident RA
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_CVD + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_CVD + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$METs, data=dat.sex))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### AD and incident RA
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_AF + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_AF + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$METs, data=dat.sex))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### CAD and incident RA
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_CAD + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_CAD + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$METs, data=dat.sex))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### HF and incident RA
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_HF + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_HF + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$METs, data=dat.sex))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### stroke and incident RA
time <-  as.numeric(dat.sex$followup_time)
status <- as.numeric(dat.sex$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_Stroke + dat.sex$Age , data = dat.sex))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.sex$exposure_Stroke + dat.sex$Age + dat.sex$BMI_group + dat.sex$Alcoholgroup + dat.sex$Smoking + dat.sex$educationgroup + dat.sex$METs, data=dat.sex))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)
rm(dat.sex)


# age group based on 60
dat.age <- subset(dat, dat$Age < 60) #260030
dat.age <- subset(dat, dat$Age >= 60)

### RA and incident CVD
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_CVD + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_CVD + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$METs, data=dat.age))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_AF + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_AF + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$METs, data=dat.age))
re2

re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident CAD
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_CAD + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_CAD + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$METs, data=dat.age))
re2

re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident HF
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_HF + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_HF + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$METs,  data=dat.age))
re2

re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident stroke
time <-  as.numeric(dat.age$followup_time)
status <- as.numeric(dat.age$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_Stroke + dat.age$sex , data = dat.age))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.age$exposure_Stroke + dat.age$sex + dat.age$BMI_group + dat.age$Alcoholgroup + dat.age$Smoking + dat.age$educationgroup + dat.age$METs, data=dat.age))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)
rm(dat.age)

# exclude the first one years
dat.ex1 <- dat[!(which(dat$followup_event == 1 & dat$followup_time <365)),] #466167

### RA and incident CVD
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_CVD + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_CVD + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$METs, data=dat.ex1))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident AF
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_AF + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_AF + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$METs, data=dat.ex1))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident CAD
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_CAD + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_CAD + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$METs, data=dat.ex1))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident HF
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_HF + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_HF + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI_group + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$METs, data=dat.ex1))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_Stroke + dat.ex1$sex + dat.ex1$Age, data = dat.ex1))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.ex1$exposure_Stroke + dat.ex1$sex + dat.ex1$Age + dat.ex1$BMI + dat.ex1$Alcoholgroup + dat.ex1$Smoking + dat.ex1$educationgroup + dat.ex1$METs, data=dat.ex1))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)

rm(dat.ex1)


# exclude kindship
dat.exkin <- subset(dat, dat$kinship == 0) #311532

### RA and incident CVD
time <-  as.numeric(dat.exkin$followup_time)
status <- as.numeric(dat.exkin$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_CVD + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_CVD + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$METs, data=dat.exkin))
re2

re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident AF
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_AF + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_AF + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$METs, data=dat.exkin))
re2
re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident CAD
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_CAD + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_CAD + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$METs, data=dat.exkin))
re2
re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident HF
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)
re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_HF + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_HF + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$METs, data=dat.exkin))
re2
re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident stroke
time <-  as.numeric(dat.ex1$followup_time)
status <- as.numeric(dat.ex1$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_Stroke + dat.exkin$sex + dat.exkin$Age, data = dat.exkin))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.exkin$exposure_Stroke + dat.exkin$sex + dat.exkin$Age + dat.exkin$BMI_group + dat.exkin$Alcoholgroup + dat.exkin$Smoking + dat.exkin$educationgroup + dat.exkin$METs, data=dat.exkin))
re2
re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)


# propensity score
library(MatchIt) 
library(survival)
dat.propensity <- dat[,c("ID","sex","Age","Townsend depreivation index","exposure_CVD","followup_event")]#441096
dat.propensity.nomissing <- dat.propensity[complete.cases(dat.propensity), ] #465734

dat.propensity.match <- matchit(exposure_CVD ~ sex + Age + `Townsend depreivation index`, data = dat.propensity.nomissing,  distance = "glm", method = "nearest", ratio = 4)
dat.propensity.match.data <- match.data(dat.propensity.match)#398230

dat.propensity.match.data <- as.data.frame(dat.propensity.match.data$ID)

dat.propensity.match.data <- merge(dat.propensity.match.data, dat, by.x = "dat.propensity.match.data$ID", by.y = "ID")

### RA and incident CVD
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_CVD + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_CVD + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$METs, data=dat.propensity.match.data))
re2
re.CVD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

### RA and incident AF
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)


re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_AF + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_AF + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$METs, data=dat.propensity.match.data))
re2

re.AF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])


### RA and incident CAD
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_CAD + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_CAD + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$METs, data=dat.propensity.match.data))
re2

re.CAD <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident HF
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_HF + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_HF + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$METs, data=dat.propensity.match.data))
re2

re.HF <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])



### RA and incident stroke
time <-  as.numeric(dat.propensity.match.data$followup_time)
status <- as.numeric(dat.propensity.match.data$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_Stroke + dat.propensity.match.data$Age + dat.propensity.match.data$sex, data = dat.propensity.match.data))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat.propensity.match.data$exposure_Stroke + dat.propensity.match.data$Age + dat.propensity.match.data$sex + dat.propensity.match.data$BMI_group + dat.propensity.match.data$Alcoholgroup + dat.propensity.match.data$Smoking + dat.propensity.match.data$educationgroup + dat.propensity.match.data$METs, data=dat.propensity.match.data))
re2

re.Stroke <- cbind(t(re1$conf.int[1,c(1,3,4)]), re1$coefficients[1,5],t(re2$conf.int[1,c(1,3,4)]), re2$coefficients[1,5])

re <- rbind(re.CVD, re.AF, re.CAD, re.HF, re.Stroke)



