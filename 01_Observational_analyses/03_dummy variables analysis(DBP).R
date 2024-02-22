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


################################################################################
# Step 1: only include European individuals
dat <- subset(dat, Enthnic=="1001"|Enthnic=="1002"|Enthnic=="1003"|Enthnic=="1")#472609

# Step 2: exclude individuals with prevalent CVD (stroke, HF, AF, CAD) at baseline
dat <- subset(dat, dat$CVD_difftime > 0)#442657
dat$CVD_selfreport <- ifelse(dat$Stroke_selfreport == 1 | dat$CAD_selfreport == 1 | dat$AF_selfreport == 1 | dat$HF_selfreport == 1, 1,0)
dat <- subset(dat, dat$CVD_selfreport == 0) #441096 some self-reported cases did not provide time, exclude them


dat$exposure_RA <- ifelse(dat$RA ==1 & dat$RA_difftime < dat$CVD_difftime , 1, 0)
table(dat$exposure_RA)

#BMI
dat$BMI_group <- ifelse(dat$BMI >30, 4,
                        ifelse(dat$BMI >25 & dat$BMI <=29.9,3,
                               ifelse(dat$BMI > 18.5 & dat$BMI <=24.9,2,1)))


# smoking
dat$Smoking[dat$Smoking== -3] <- 9

#alcohol consumption#Ref1:Never
dat <- dat %>% mutate(Alcoholgroup=if_else(Alcohol=="6",1,
                                           ifelse(Alcohol=="5",2,
                                                  ifelse(Alcohol=="4",3,
                                                         ifelse(Alcohol=="3",4,
                                                                ifelse(Alcohol=="2",5,
                                                                       ifelse(Alcohol=="1",6, 9)))))))
#Education
dat <- dat %>% mutate(educationgroup=ifelse(Education.qua == "1" | Education.qua == "6",3,
                                            ifelse(Education.qua == "2"| Education.qua == "3"| Education.qua == "4"|Education.qua == "5",2,
                                                   ifelse(Education.qua == "-7",1, 9))))
# hypertension
dat$Hyp_baseline <- ifelse(dat$Hyp == 1 & dat$Hyp_difftime < 0 , 1,0) # baseline Hypertension

#T2D
dat$T2D_baseline <- ifelse(dat$T2D == 1 & dat$T2D_difftime < 0 , 1,0) # baseline Hypertension


# DBP
dat$DBP.group <- ifelse(dat$diastolic.blood.pressure >= 90,1,0)
dat <- dat %>% mutate(RADBP =ifelse(exposure_RA == 0 & DBP.group == 0,0,
                                            ifelse(exposure_RA == 1 & DBP.group == 0,2,
                                                   ifelse(exposure_RA == 0 & DBP.group == 1,1,3))))




################################################################################
# Number of event and individuals
table(dat$RADBP, dat$followup_event)

# generate the plot
time <- dat$followup_time/365
library(RColorBrewer)
library(ggpubr)
library(survminer)
library(survival)
colors <- brewer.pal(8,"Dark2")[c(1:4)]

ggsurv <- ggsurvplot(survfit(Surv(time,dat$followup_event) ~ dat$RADBP, data = dat),
                     pval = F,
                     palette = colors,
                     legend.title=" ",
                     legend.labs = c("RA(-) High DBP(-)", "RA(+) High DBP(-)","RA(-) High DBP(+)","RA(+) High DBP(+)"),
                     legend = c(0.25,0.8),#legend 
                     fun = "cumhaz",
                     xlab = "Follow-up(years)",
                     ylim = c(0,0.4)
)


ggsurv$plot <- ggsurv$plot +
  theme(legend.title = element_text(size = 2),legend.text =  element_text(size=10, color ="black",face="bold"), axis.text.x = element_text(size= 12, color="black" ,face="bold"),axis.text.y = element_text(size= 12, color="black" ,face="bold")) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsurv

### RA and incident CVD
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ as.factor(dat$RADBP) + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ as.factor(dat$RADBP) + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup + dat$Hyp_baseline + dat$T2D_baseline + dat$METs + dat$TC, data=dat))
re2

#################################################################################
#CVD-RA
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

dat$exposure_CVD <- ifelse(dat$CVD == 1 & dat$CVD_difftime < dat$RA_difftime , 1, 0)

#################################
#BMI
dat$BMI_group <- ifelse(dat$BMI >30, 4,
                        ifelse(dat$BMI >25 & dat$BMI <=29.9,3,
                               ifelse(dat$BMI > 18.5 & dat$BMI <=24.9,2,1)))


# smoking
dat$Smoking[dat$Smoking== -3] <- 9

#alcohol consumption#Ref1:Never
dat <- dat %>% mutate(Alcoholgroup=if_else(Alcohol=="6",1,
                                           ifelse(Alcohol=="5",2,
                                                  ifelse(Alcohol=="4",3,
                                                         ifelse(Alcohol=="3",4,
                                                                ifelse(Alcohol=="2",5,
                                                                       ifelse(Alcohol=="1",6, 9)))))))
#Education
dat <- dat %>% mutate(educationgroup=ifelse(Education.qua == "1" | Education.qua == "6",3,
                                            ifelse(Education.qua == "2"| Education.qua == "3"| Education.qua == "4"|Education.qua == "5",2,
                                                   ifelse(Education.qua == "-7",1, 9))))
# hypertension
dat$Hyp_baseline <- ifelse(dat$Hyp == 1 & dat$Hyp_difftime < 0 , 1,0) # baseline Hypertension

#T2D
dat$T2D_baseline <- ifelse(dat$T2D == 1 & dat$T2D_difftime < 0 , 1,0) # baseline Hypertension


# DBP
dat$DBP.group <- ifelse(dat$diastolic.blood.pressure >= 90,1,0)
dat <- dat %>% mutate(CVDDBP =ifelse(exposure_CVD == 0 & DBP.group == 0,0,
                                    ifelse(exposure_CVD == 1 & DBP.group == 0,2,
                                           ifelse(exposure_CVD == 0 & DBP.group == 1,1,3))))



################################################################################
# Number of event and individuals
table(dat$CVDDBP, dat$followup_event)

# generate the plot
time <- dat$followup_time/365
library(RColorBrewer)
library(ggpubr)
library(survminer)
library(survival)
colors <- brewer.pal(8,"Dark2")[c(1:4)]

ggsurv <- ggsurvplot(survfit(Surv(time,dat$followup_event) ~ dat$CVDDBP, data = dat),
                     pval = F,
                     palette = colors,
                     legend.title=" ",
                     legend.labs = c("CVD(-) High DBP(-)", "CVD(+) High DBP(-)","CVD(-) High DBP(+)","CVD(+) High DBP(+)"),
                     legend = c(0.25,0.8),#legend 
                     fun = "cumhaz",
                     xlab = "Follow-up(years)",
                     ylim = c(0,0.03)
)


ggsurv$plot <- ggsurv$plot +
  theme(legend.title = element_text(size = 2),legend.text =  element_text(size=10, color ="black",face="bold"), axis.text.x = element_text(size= 12, color="black" ,face="bold"),axis.text.y = element_text(size= 12, color="black" ,face="bold")) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

ggsurv


### CVD and incident RA
time <-  as.numeric(dat$followup_time)
status <- as.numeric(dat$followup_event)

re1 <- summary(coxph(Surv(time, status) ~ as.factor(dat$CVDDBP) + dat$Age + dat$sex , data = dat))
re1
re2 <- summary(coxph(Surv(time, status) ~ dat$CVDDBP + dat$Age + dat$sex + dat$BMI_group + dat$Alcoholgroup + dat$Smoking + dat$educationgroup  + dat$METs , data=dat))
re2





