library(data.table)
library(tidyverse)
library(stringr)
library(lubridate)

setwd("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB")
dat <- fread("UKBB_pbc.tab")
dat.sup <- fread("UKBB_pbc_supp.tab")
dat.sup2 <- fread("UKBB_pbc_supp2.tab")
dat.sup3 <- fread("UKBB_pbc_supp3.tab")#death time and death cause

dat <- merge(dat, dat.sup, by = "f.eid")
dat <- merge(dat, dat.sup2, by = "f.eid")
dat <- merge(dat, dat.sup3, by = "f.eid")

rm(dat.sup)
rm(dat.sup2)
rm(dat.sup3)

dat <- as.data.frame(dat)


#################################################################
#loss time: 191
dat.followup <- dat[,c("f.eid","f.191.0.0")]
colnames(dat.followup)[2] <- "losetime"

# enter time: 53
dat.entertime <- dat[,c("f.eid","f.53.0.0","f.53.1.0","f.53.2.0","f.53.3.0")]
entertime <- as.Date(vector())
for(i in 1:nrow(dat.entertime)){
  
  if(complete.cases(dat.entertime$f.53.0.0[i]))
  {
    entertime[i]=dat.entertime$f.53.0.0[i]
    
  }
  else if (complete.cases(dat.entertime$f.53.1.0[i]))
  {
    entertime[i]=dat.entertime$f.53.1.0[i]
    
  }
  else if (complete.cases(dat.entertime$f.53.2.0[i]))
  {
    entertime[i]=dat.entertime$f.53.2.0[i]
  }
  else if (complete.cases(dat.entertime$f.53.3.0[i]))
  {
    entertime[i]=dat.entertime$f.53.3.0[i]
  }
  
  else {
    entertime[i]=dat.entertime$f.53.0.0[i]
    
  }
} 

dat.entertime <- cbind(dat.entertime, entertime)
dat.followup <- merge(dat.followup, dat.entertime[,c("f.eid","entertime")], by.x = "f.eid", by.y = "f.eid")
rm(dat.entertime)

## end time
dat.residence <- dat[,c(1,392)] #Field 20118
dat.residence$redince <- ifelse(dat.residence$f.20118.0.0==1| dat.residence$f.20118.0.0==2 |dat.residence$f.20118.0.0==3 | dat.residence$f.20118.0.0==4 |dat.residence$f.20118.0.0==5|dat.residence$f.20118.0.0==6 |dat.residence$f.20118.0.0==7 |dat.residence$f.20118.0.0==8,"England",
                                ifelse(dat.residence$f.20118.0.0==11 |dat.residence$f.20118.0.0==12 | dat.residence$f.20118.0.0==13 |dat.residence$f.20118.0.0==14 |dat.residence$f.20118.0.0==15|dat.residence$f.20118.0.0==16|dat.residence$f.20118.0.0==17|dat.residence$f.20118.0.0==18,"Scotland", NA))
dat.followup <- merge(dat.followup , dat.residence[,c(1,3)], by.x="f.eid", by.y="f.eid")

dat.followup$endtime <- ifelse(dat.followup$redince == "England", as.Date("2020-02-29"),
                               ifelse(dat.followup$redince =="Scotland", as.Date("2021-01-31"), 0))
dat.followup$endtime <- as.Date(dat.followup$endtime, "%Y-%m-%d",origin="1970-01-01")


# death time and death cause (Underlying (primary) cause of death: ICD10)
dat.death <- dat[,c("f.eid","f.40001.0.0","f.40001.1.0")]
dat.death.time <- dat[,c("f.eid","f.40000.0.0","f.40000.1.0")]

## extract the major outcome (cardiovascular disease and RA)
### RA
RA_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/RA_ICD10.csv",head=F)
v <- as.vector(RA_ICD10_code$V1)
#### extract the location which contain ICD10 code of RA
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.death == x,arr.ind = T)))
Loc_row <- rowcol[,1]
rowcol <- as.data.frame(cbind(dat.death[Loc_row,1],rowcol))

RA_ICD10 <- dat.death[Loc_row,] #45
RA_ICD10 <- RA_ICD10[!duplicated(RA_ICD10$f.eid),] #45
RA_ICD10_time <- as.data.frame(dat.death.time[Loc_row,])
RA_ICD10_time <- RA_ICD10_time[!duplicated(RA_ICD10_time$f.eid),]#830

for(i in 1:nrow(RA_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == RA_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.death.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  RA_ICD10_time[i,4] <- sort.temp[1]
}

dat.followup$RA_death <- ifelse(dat.followup$f.eid %in% RA_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, RA_ICD10_time[,c(1,4)], by.x = "f.eid", by.y = "f.eid", all=T)
colnames(dat.followup)[5] <- "RA_death_time"

### Stroke
Stroke_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/stroke_ICD10.csv",head=F)
v <- as.vector(Stroke_ICD10_code$V1)
#### extract the location which contain ICD10 code of stroke
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.death == x,arr.ind = T)))
Loc_row <- rowcol[,1]
rowcol <- as.data.frame(cbind(dat.death[Loc_row,1],rowcol))

Stroke_ICD10 <- dat.death[Loc_row,] #1204
Stroke_ICD10 <- Stroke_ICD10[!duplicated(Stroke_ICD10$f.eid),] #1204
Stroke_ICD10_time <- as.data.frame(dat.death.time[Loc_row,])
Stroke_ICD10_time <- Stroke_ICD10_time[!duplicated(Stroke_ICD10_time$f.eid),]#1204

for(i in 1:nrow(Stroke_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == Stroke_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.death.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  Stroke_ICD10_time[i,4] <- sort.temp[1]
}

dat.followup$Stroke_death <- ifelse(dat.followup$f.eid %in% Stroke_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, Stroke_ICD10_time[,c(1,4)], by.x = "f.eid", by.y = "f.eid", all=T)
colnames(dat.followup)[7] <- "Stroke_death_time"


### CAD
CAD_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/CAD_ICD10.csv",head=F)
v <- as.vector(CAD_ICD10_code$V1)
#### extract the location which contain ICD10 code of CAD
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.death == x,arr.ind = T)))
Loc_row <- rowcol[,1]
rowcol <- as.data.frame(cbind(dat.death[Loc_row,1],rowcol))

CAD_ICD10 <- dat.death[Loc_row,] #3718
CAD_ICD10 <- CAD_ICD10[!duplicated(CAD_ICD10$f.eid),] #3712
CAD_ICD10_time <- as.data.frame(dat.death.time[Loc_row,])#3718
CAD_ICD10_time <- CAD_ICD10_time[!duplicated(CAD_ICD10_time$f.eid),]#3712

for(i in 1:nrow(CAD_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == CAD_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.death.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  CAD_ICD10_time[i,4] <- sort.temp[1]
}

dat.followup$CAD_death <- ifelse(dat.followup$f.eid %in% CAD_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, CAD_ICD10_time[,c(1,4)], by.x = "f.eid", by.y = "f.eid", all=T)
colnames(dat.followup)[9] <- "CAD_death_time"

### AF
AF_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/AF_ICD10.csv",head=F)
v <- as.vector(AF_ICD10_code$V1)
#### extract the location which contain ICD10 code of AF
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.death == x,arr.ind = T)))
Loc_row <- rowcol[,1]
rowcol <- as.data.frame(cbind(dat.death[Loc_row,1],rowcol))

AF_ICD10 <- dat.death[Loc_row,] #131
AF_ICD10 <- AF_ICD10[!duplicated(AF_ICD10$f.eid),] #131
AF_ICD10_time <- as.data.frame(dat.death.time[Loc_row,])#131
AF_ICD10_time <- AF_ICD10_time[!duplicated(AF_ICD10_time$f.eid),]#131

for(i in 1:nrow(AF_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == AF_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.death.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  AF_ICD10_time[i,4] <- sort.temp[1]
}

dat.followup$AF_death <- ifelse(dat.followup$f.eid %in% AF_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, AF_ICD10_time[,c(1,4)], by.x = "f.eid", by.y = "f.eid", all=T)
colnames(dat.followup)[11] <- "AF_death_time"


### HF
HF_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/HF_ICD10.csv",head=F)
v <- as.vector(HF_ICD10_code$V1)
#### extract the location which contain ICD10 code of AF
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.death == x,arr.ind = T)))
Loc_row <- rowcol[,1]
rowcol <- as.data.frame(cbind(dat.death[Loc_row,1],rowcol))

HF_ICD10 <- dat.death[Loc_row,] #128
HF_ICD10 <- HF_ICD10[!duplicated(HF_ICD10$f.eid),] #127
HF_ICD10_time <- as.data.frame(dat.death.time[Loc_row,])#128
HF_ICD10_time <- HF_ICD10_time[!duplicated(HF_ICD10_time$f.eid),]#131

for(i in 1:nrow(HF_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == HF_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.death.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  HF_ICD10_time[i,4] <- sort.temp[1]
}

dat.followup$HF_death <- ifelse(dat.followup$f.eid %in% HF_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, HF_ICD10_time[,c(1,4)], by.x = "f.eid", by.y = "f.eid", all=T)
colnames(dat.followup)[13] <- "HF_death_time"

####################################### ICD10 ########################################
## ICD 10 code: 41270
dat.icd10 <- dat[,c(1,479:704)]
## ICD 10 time: 41280
dat.icd10.time <- dat[,c(1,752:977)]

###RA
RA_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/RA_ICD10.csv",head=F)
v <- as.vector(RA_ICD10_code$V1)
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

RA_ICD10 <- dat.icd10[Loc_row,]#16367
RA_ICD10 <- RA_ICD10[!duplicated(RA_ICD10$f.eid),]#9163

RA_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])
RA_ICD10_time <- RA_ICD10_time[!duplicated(RA_ICD10_time$f.eid),]

for(i in 1:nrow(RA_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == RA_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  RA_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$RA_ICD10 <- ifelse(dat.followup$f.eid %in% RA_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, RA_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[15] <- "RA_ICD10_time"



###stroke
Stroke_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/stroke_ICD10.csv",head=F)
v <- as.vector(Stroke_ICD10_code$V1)
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

Stroke_ICD10 <- dat.icd10[Loc_row,] #15021
Stroke_ICD10 <- Stroke_ICD10[!duplicated(Stroke_ICD10$f.eid),] #12387

Stroke_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])#15021
Stroke_ICD10_time <- Stroke_ICD10_time[!duplicated(Stroke_ICD10_time$f.eid),]#12387

for(i in 1:nrow(Stroke_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == Stroke_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  Stroke_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$Stroke_ICD10 <- ifelse(dat.followup$f.eid %in% Stroke_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, Stroke_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[17] <- "Stroke_ICD10_time"

###CAD
CAD_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/CAD_ICD10.csv",head=F)
v <- as.vector(CAD_ICD10_code$V1)
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

CAD_ICD10 <- dat.icd10[Loc_row,] #112025
CAD_ICD10 <- CAD_ICD10[!duplicated(CAD_ICD10$f.eid),] #49589

CAD_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])#112025
CAD_ICD10_time <- CAD_ICD10_time[!duplicated(CAD_ICD10_time$f.eid),]#49589

for(i in 1:nrow(CAD_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == CAD_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  CAD_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$CAD_ICD10 <- ifelse(dat.followup$f.eid %in% CAD_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, CAD_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[19] <- "CAD_ICD10_time"


### AF
AF_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/AF_ICD10.csv",head=F)
v <- as.vector(AF_ICD10_code$V1)
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

AF_ICD10 <- dat.icd10[Loc_row,] #50485
AF_ICD10 <- AF_ICD10[!duplicated(AF_ICD10$f.eid),] #34642

AF_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])#50485
AF_ICD10_time <- AF_ICD10_time[!duplicated(AF_ICD10_time$f.eid),]#49589

for(i in 1:nrow(AF_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == AF_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  AF_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$AF_ICD10 <- ifelse(dat.followup$f.eid %in% AF_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, AF_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[21] <- "AF_ICD10_time"


### HF
HF_ICD10_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/HF_ICD10.csv",head=F)
v <- as.vector(HF_ICD10_code$V1)
rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

HF_ICD10 <- dat.icd10[Loc_row,] #21536
HF_ICD10 <- HF_ICD10[!duplicated(HF_ICD10$f.eid),] #15867

HF_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])#21536
HF_ICD10_time <- HF_ICD10_time[!duplicated(HF_ICD10_time$f.eid),]#15867

for(i in 1:nrow(HF_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == HF_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  HF_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$HF_ICD10 <- ifelse(dat.followup$f.eid %in% HF_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, HF_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[23] <- "HF_ICD10_time"

### Hypertension
v = as.vector(c("I10"))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

Hyp_ICD10 <- dat.icd10[Loc_row,] #146954
Hyp_ICD10 <- Hyp_ICD10[!duplicated(Hyp_ICD10$f.eid),] #146954

Hyp_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])#146954
Hyp_ICD10_time <- Hyp_ICD10_time[!duplicated(Hyp_ICD10_time$f.eid),]#146954

for(i in 1:nrow(Hyp_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == Hyp_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  Hyp_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$Hyp_ICD10 <- ifelse(dat.followup$f.eid %in% Hyp_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, Hyp_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[25] <- "Hyp_ICD10_time"


### Type 2 diabetes
v <- as.vector(unlist(c("E11","E110","E111","E112","E113","E114","E115","E116","E117","E118","E119")))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd10 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd10[Loc_row,1],rowcol))

T2D_ICD10 <- dat.icd10[Loc_row,] 
T2D_ICD10 <- T2D_ICD10[!duplicated(T2D_ICD10$f.eid),] 

T2D_ICD10_time <- as.data.frame(dat.icd10.time[Loc_row,])
T2D_ICD10_time <- T2D_ICD10_time[!duplicated(T2D_ICD10_time$f.eid),]#39564

for(i in 1:nrow(T2D_ICD10_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == T2D_ICD10_time$f.eid[i])
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.icd10.time[e,c]
    temp.all[j,1] <- format(as.Date(temp))
  }
  sort.temp <- sort(temp.all$NA.)
  T2D_ICD10_time[i,228] <- sort.temp[1]
}

dat.followup$T2D_ICD10 <- ifelse(dat.followup$f.eid %in% T2D_ICD10$f.eid,1,0)
dat.followup <- merge(dat.followup, T2D_ICD10_time[,c(1,228)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[27] <- "T2D_ICD10_time"



rm(dat.icd10)
rm(dat.icd10.time)
###############################################################################
###################### ICD9 ##################################################
## ICD9:40013
dat.icd9 <- dat[,c(1,442:456)]
## ICD9 time: 41281
dat.icd9.time <- as.data.frame(dat[,c(1,978:1024)])

### RA
RA_ICD9_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/RA_ICD9.csv",head=F)
v <- as.vector(unlist(RA_ICD9_code[,1]))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

RA_ICD9 <- dat.icd9[Loc_row,]#0
RA_ICD9 <- RA_ICD9[!duplicated(RA_ICD9$f.eid),]#0

RA_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
RA_ICD9_time <- RA_ICD9_time[!duplicated(RA_ICD9_time$f.eid),]#0


dat.followup$RA_ICD9 <- 0
dat.followup$RA_ICD9_time <- NA

### stroke
Stroke_ICD9_code <- read.csv("C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/disease code/stroke_ICD9.csv",head=F)
v <- as.vector(Stroke_ICD9_code$V1)

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

Stroke_ICD9 <- dat.icd9[Loc_row,]#0
Stroke_ICD9 <- Stroke_ICD9[!duplicated(Stroke_ICD9$f.eid),]#0

Stroke_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
Stroke_ICD9_time <- Stroke_ICD9_time[!duplicated(Stroke_ICD9_time$f.eid),]#0

dat.followup$Stroke_ICD9 <- 0
dat.followup$Stroke_ICD9_time <- NA

### CAD
CAD_ICD9_code <- as.vector(unlist(c("410","4109","411","4119","412","4129","4140","4148","4149")))
v <- as.vector(CAD_ICD9_code)

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

CAD_ICD9 <- dat.icd9[Loc_row,]#0
CAD_ICD9 <- CAD_ICD9[!duplicated(CAD_ICD9$f.eid),]#0

CAD_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
CAD_ICD9_time <- CAD_ICD9_time[!duplicated(CAD_ICD9_time$f.eid),]#0

dat.followup$CAD_ICD9 <- 0
dat.followup$CAD_ICD9_time <- NA

### AF
v = as.vector(unlist(c("4273")))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

AF_ICD9 <- dat.icd9[Loc_row,]#0
AF_ICD9 <- AF_ICD9[!duplicated(AF_ICD9$f.eid),]#0

AF_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
AF_ICD9_time <- AF_ICD9_time[!duplicated(AF_ICD9_time$f.eid),]#0

dat.followup$AF_ICD9 <- 0
dat.followup$AF_ICD9_time <- NA

### HF
v = as.vector(unlist(c("428","4280","4281")))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

HF_ICD9 <- dat.icd9[Loc_row,]#0
HF_ICD9 <- HF_ICD9[!duplicated(HF_ICD9$f.eid),]#0

HF_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
HF_ICD9_time <- HF_ICD9_time[!duplicated(HF_ICD9_time$f.eid),]#0

dat.followup$HF_ICD9 <- 0
dat.followup$HF_ICD9_time <- NA




### hypertension
v = as.vector(unlist(c("401","4010","4011","4019")))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

Hyp_ICD9 <- dat.icd9[Loc_row,]#0
Hyp_ICD9 <- Hyp_ICD9[!duplicated(Hyp_ICD9$f.eid),]#0

Hyp_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
Hyp_ICD9_time <- Hyp_ICD9_time[!duplicated(Hyp_ICD9_time$f.eid),]#0

dat.followup$Hyp_ICD9 <- 0
dat.followup$Hyp_ICD9_time <- NA


### T2D
v = as.vector(unlist(c("2500","25000","25001","25009","2501","25010","25011","25019","2502","25029","2503","2504","2505","2509","25099")))

rowcol <- do.call('rbind',lapply(v,function(x) which(dat.icd9 == x,arr.ind = T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.icd9[Loc_row,1],rowcol))

T2D_ICD9 <- dat.icd9[Loc_row,]#0
T2D_ICD9 <- T2D_ICD9[!duplicated(T2D_ICD9$f.eid),]#0

T2D_ICD9_time <- as.data.frame(dat.icd9.time[Loc_row,])#0
T2D_ICD9_time <- T2D_ICD9_time[!duplicated(T2D_ICD9_time$f.eid),]#0

dat.followup$T2D_ICD9 <- 0
dat.followup$T2D_ICD9_time <- NA


rm(dat.icd9)
rm(dat.icd9.time)
################################################################################
#################### self-reported #########################################
## self-reported code: 20002
dat.selfreport <- dat[,c(1,60:168)]
## self-reported time:20008
dat.selfreport.time <- as.data.frame(dat[,c(1,1076:1211)])

### RA
v = as.vector(unlist(c("1464")))

rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))

RA_selfreport <- dat.selfreport[Loc_row,]#6301
RA_selfreport <- RA_selfreport[!duplicated(RA_selfreport$f.eid),]#5985

RA_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#6301
RA_selfreport_time <- RA_selfreport_time[!duplicated(RA_selfreport_time$f.eid),]#5985

for(i in 1:nrow(RA_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == RA_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  RA_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$RA_selfreport <- ifelse(dat.followup$f.eid %in% RA_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, RA_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[43] <- "RA_selfreport_time"

### Stroke
v = as.vector(unlist(c("1081","1086","1491","1583")))
rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))

Stroke_selfreport <- dat.selfreport[Loc_row,]#8217
Stroke_selfreport <- Stroke_selfreport[!duplicated(Stroke_selfreport$f.eid),]#7710

Stroke_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#8217
Stroke_selfreport_time <- Stroke_selfreport_time[!duplicated(Stroke_selfreport_time$f.eid),]#7710

for(i in 1:nrow(Stroke_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == Stroke_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  Stroke_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$Stroke_selfreport <- ifelse(dat.followup$f.eid %in% Stroke_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, Stroke_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[45] <- "Stroke_selfreport_time"

### CAD
v = as.vector(unlist(c("1075")))

rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))

CAD_selfreport <- dat.selfreport[Loc_row,]#12997
CAD_selfreport <- CAD_selfreport[!duplicated(CAD_selfreport$f.eid),]#12091

CAD_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#12997
CAD_selfreport_time <- CAD_selfreport_time[!duplicated(CAD_selfreport_time$f.eid),]#12091

for(i in 1:nrow(CAD_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == CAD_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  CAD_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$CAD_selfreport <- ifelse(dat.followup$f.eid %in% CAD_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, CAD_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[47] <- "CAD_selfreport_time"


### AF
v = as.vector(unlist(c("1471")))
rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))

AF_selfreport <- dat.selfreport[Loc_row,]#4831
AF_selfreport <- AF_selfreport[!duplicated(AF_selfreport$f.eid),]#4474

AF_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#4831
AF_selfreport_time <- AF_selfreport_time[!duplicated(AF_selfreport_time$f.eid),]#4474

for(i in 1:nrow(AF_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == AF_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  AF_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$AF_selfreport <- ifelse(dat.followup$f.eid %in% AF_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, AF_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[49] <- "AF_selfreport_time"


### HF
v = as.vector(unlist(c("1076")))
rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))

HF_selfreport <- dat.selfreport[Loc_row,]#380
HF_selfreport <- HF_selfreport[!duplicated(HF_selfreport$f.eid),]#370

HF_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#380
HF_selfreport_time <- HF_selfreport_time[!duplicated(HF_selfreport_time$f.eid),]#370

for(i in 1:nrow(HF_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == HF_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  HF_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$HF_selfreport <- ifelse(dat.followup$f.eid %in% HF_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, HF_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[51] <- "HF_selfreport_time"

### hypertension
v = as.vector(unlist(c("1065")))
rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))


Hyp_selfreport <- dat.selfreport[Loc_row,]#147202
Hyp_selfreport <- Hyp_selfreport[!duplicated(Hyp_selfreport$f.eid),]#136041

Hyp_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#147202
Hyp_selfreport_time <- Hyp_selfreport_time[!duplicated(Hyp_selfreport_time$f.eid),]#136041

for(i in 1:nrow(Hyp_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == Hyp_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  Hyp_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$Hyp_selfreport <- ifelse(dat.followup$f.eid %in% Hyp_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, Hyp_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[53] <- "Hyp_selfreport_time"

### T2D
v = as.vector(unlist(c("1223")))
rowcol<-do.call('rbind',lapply(v,function(x) which(dat.selfreport == x,arr.ind=T)))
Loc_row<-rowcol[,1]
rowcol <- as.data.frame(cbind(dat.selfreport[Loc_row,1],rowcol))


T2D_selfreport <- dat.selfreport[Loc_row,]#4622
T2D_selfreport <- T2D_selfreport[!duplicated(T2D_selfreport$f.eid),]#4485

T2D_selfreport_time <- as.data.frame(dat.selfreport.time[Loc_row,])#4622
T2D_selfreport_time <- T2D_selfreport_time[!duplicated(T2D_selfreport_time$f.eid),]#4485

for(i in 1:nrow(T2D_selfreport_time)){
  temp <- data.frame(NA)
  temp.all <- data.frame(NA)
  loc <- subset(rowcol,rowcol$V1 == T2D_selfreport_time$f.eid[i]) #
  a_row <- as.data.frame(loc[,2])
  b_col <- as.data.frame(loc[,3])
  for (j in 1:nrow(b_col)){
    e <- as.numeric(a_row[j,])
    c <- as.numeric(b_col[j,])
    temp <- dat.selfreport.time[e,c]#
    temp.all[j,1] <- temp
  }
  sort.temp <- sort(temp.all$NA.)
  T2D_selfreport_time[i,138] <- sort.temp[1]#
}

dat.followup$T2D_selfreport <- ifelse(dat.followup$f.eid %in% T2D_selfreport$f.eid,1,0)
dat.followup <- merge(dat.followup, T2D_selfreport_time[,c(1,138)], by.x = "f.eid", by.y ="f.eid",all = T)
colnames(dat.followup)[55] <- "T2D_selfreport_time"

write.csv(dat.followup, "ukbb_RA&CAD_followup_1.csv",row.names = F)

################################################################
################################################################
dat.followup <- fread("ukbb_RA&CAD_followup_1.csv")


#ICD 10
dat.followup$RA_ICD10_difftime <- as.numeric(difftime(dat.followup$RA_ICD10_time , dat.followup$entertime,units="days"), units="days")
dat.followup$Stroke_ICD10_difftime <- as.numeric(difftime(dat.followup$Stroke_ICD10_time , dat.followup$entertime,units="days"), units="days")
dat.followup$CAD_ICD10_difftime <- as.numeric(difftime(dat.followup$CAD_ICD10_time , dat.followup$entertime, units="days"),units="days")
dat.followup$AF_ICD10_difftime <- as.numeric(difftime(dat.followup$AF_ICD10_time , dat.followup$entertime,units="days"), units="days")
dat.followup$HF_ICD10_difftime <- as.numeric(difftime(dat.followup$HF_ICD10_time , dat.followup$entertime, units="days"),units="days")
dat.followup$Hyp_ICD10_difftime <- as.numeric(difftime(dat.followup$Hyp_ICD10_time , dat.followup$entertime, units="days"),units="days")
dat.followup$T2D_ICD10_difftime <- as.numeric(difftime(dat.followup$T2D_ICD10_time , dat.followup$entertime, units="days"),units="days")


#ICD 9
dat.followup$RA_ICD9_difftime <- as.numeric(difftime(dat.followup$RA_ICD9_time , dat.followup$entertime,units="days"), units="days")
dat.followup$Stroke_ICD9_difftime <- as.numeric(difftime(dat.followup$Stroke_ICD9_time , dat.followup$entertime,units="days"), units="days")
dat.followup$CAD_ICD9_difftime <- as.numeric(difftime(dat.followup$CAD_ICD9_time , dat.followup$entertime, units="days"),units="days")
dat.followup$AF_ICD9_difftime <- as.numeric(difftime(dat.followup$AF_ICD9_time , dat.followup$entertime,units="days"), units="days")
dat.followup$HF_ICD9_difftime <- as.numeric(difftime(dat.followup$HF_ICD9_time , dat.followup$entertime, units="days"),units="days")
dat.followup$Hyp_ICD9_difftime <- as.numeric(difftime(dat.followup$Hyp_ICD9_time , dat.followup$entertime,units="days"), units="days")
dat.followup$T2D_ICD9_difftime <- as.numeric(difftime(dat.followup$T2D_ICD9_time , dat.followup$entertime,units="days"), units="days")

#self-reported

dat.followup$a <- trunc(dat.followup$RA_selfreport_time)##
dat.followup$b <- trunc((dat.followup$RA_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$RA_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$RA_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")

dat.followup$a <- trunc(dat.followup$Stroke_selfreport_time)##
dat.followup$b <- trunc((dat.followup$Stroke_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$Stroke_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$Stroke_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


dat.followup$a <- trunc(dat.followup$CAD_selfreport_time)##
dat.followup$b <- trunc((dat.followup$CAD_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$CAD_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$CAD_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


dat.followup$a <- trunc(dat.followup$AF_selfreport_time)##
dat.followup$b <- trunc((dat.followup$AF_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$AF_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$AF_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


dat.followup$a <- trunc(dat.followup$HF_selfreport_time)##
dat.followup$b <- trunc((dat.followup$HF_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$HF_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$HF_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


dat.followup$a <- trunc(dat.followup$Hyp_selfreport_time)##
dat.followup$b <- trunc((dat.followup$Hyp_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$Hyp_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$Hyp_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


dat.followup$a <- trunc(dat.followup$T2D_selfreport_time)##
dat.followup$b <- trunc((dat.followup$T2D_selfreport_time - dat.followup$a)*12)##
dat.followup$c <- as.numeric(dat.followup$T2D_selfreport_time) - as.numeric(dat.followup$a) - as.numeric(dat.followup$b/12) ##
dat.followup$d <- ceiling(dat.followup$c*365)
dat.followup$b[dat.followup$b == 0] <- 1
dat.followup$d[dat.followup$d == 0] <- 1
dat.followup$e <- as.Date(paste(dat.followup$a , dat.followup$b , dat.followup$d),"%Y%m%d")
dat.followup$T2D_selfreport_difftime <- as.numeric(difftime(dat.followup$e , dat.followup$entertime,units="days"), units="days")


# death time
dat.followup$RA_death_difftime <- as.numeric(difftime(dat.followup$RA_death_time , dat.followup$entertime,units="days"), units="days")
dat.followup$Stroke_death_difftime <- as.numeric(difftime(dat.followup$Stroke_death_time , dat.followup$entertime,units="days"), units="days")
dat.followup$CAD_death_difftime <- as.numeric(difftime(dat.followup$CAD_death_time , dat.followup$entertime, units="days"),units="days")
dat.followup$AF_death_difftime <- as.numeric(difftime(dat.followup$AF_death_time , dat.followup$entertime,units="days"), units="days")
dat.followup$HF_death_difftime <- as.numeric(difftime(dat.followup$HF_death_time , dat.followup$entertime, units="days"),units="days")


# extract the minimum of time across ICD 10, ICD 9 and death time as the outcome incidence time
dat.followup$RA_difftime <- apply(dat.followup[, c("RA_ICD10_difftime","RA_ICD9_difftime","RA_selfreport_difftime")],1,min, na.rm=T)
dat.followup$Stroke_difftime <- apply(dat.followup[, c("Stroke_ICD10_difftime","Stroke_ICD9_difftime","Stroke_selfreport_difftime")],1,min, na.rm=T)
dat.followup$CAD_difftime <- apply(dat.followup[, c("CAD_ICD10_difftime","CAD_ICD9_difftime","CAD_selfreport_difftime")],1,min, na.rm=T)
dat.followup$AF_difftime <- apply(dat.followup[, c("AF_ICD10_difftime","AF_ICD9_difftime","AF_selfreport_difftime")],1,min, na.rm=T)
dat.followup$HF_difftime <- apply(dat.followup[, c("HF_ICD10_difftime","HF_ICD9_difftime","HF_selfreport_difftime")],1,min, na.rm=T)
dat.followup$Hyp_difftime <- apply(dat.followup[, c("Hyp_ICD10_difftime","Hyp_ICD9_difftime","Hyp_selfreport_difftime")],1,min, na.rm=T)
dat.followup$T2D_difftime <- apply(dat.followup[, c("T2D_ICD10_difftime","T2D_ICD9_difftime","T2D_selfreport_difftime")],1,min, na.rm=T)

dat.followup$RA <- ifelse(dat.followup$RA_ICD10 == 1 | dat.followup$RA_ICD9 == 1 | dat.followup$RA_selfreport == 1,1,0)
dat.followup$Stroke <- ifelse(dat.followup$Stroke_ICD10 == 1 | dat.followup$Stroke_ICD9 ==1 | dat.followup$Stroke_selfreport ==1,1,0)
dat.followup$CAD <- ifelse(dat.followup$CAD_ICD10==1 | dat.followup$CAD_ICD9 == 1 | dat.followup$CAD_selfreport ==1 ,1,0)
dat.followup$AF <- ifelse(dat.followup$AF_ICD10==1 | dat.followup$AF_ICD9 ==1 | dat.followup$AF_selfreport ==1 ,1,0)
dat.followup$HF <- ifelse(dat.followup$HF_ICD10==1 | dat.followup$HF_ICD9 ==1 | dat.followup$HF_selfreport ==1 ,1,0)
dat.followup$Hyp <- ifelse(dat.followup$Hyp_ICD10==1 | dat.followup$Hyp_ICD9 ==1 | dat.followup$Hyp_selfreport ==1 ,1,0)
dat.followup$T2D <- ifelse(dat.followup$T2D_ICD10==1 | dat.followup$T2D_ICD9 ==1 | dat.followup$T2D_selfreport ==1,1,0)

## lose time
dat.followup$losetime_difftime <- as.numeric(difftime(dat.followup$losetime , dat.followup$entertime,units="days"), units="days")

## end time
dat.residence <- dat[,c(1,392)] #Field 20118
dat.residence$redince <- ifelse(dat.residence$f.20118.0.0==1| dat.residence$f.20118.0.0==2 |dat.residence$f.20118.0.0==3 | dat.residence$f.20118.0.0==4 |dat.residence$f.20118.0.0==5|dat.residence$f.20118.0.0==6 |dat.residence$f.20118.0.0==7 |dat.residence$f.20118.0.0==8,"England",
                                ifelse(dat.residence$f.20118.0.0==11 |dat.residence$f.20118.0.0==12 | dat.residence$f.20118.0.0==13 |dat.residence$f.20118.0.0==14 |dat.residence$f.20118.0.0==15|dat.residence$f.20118.0.0==16|dat.residence$f.20118.0.0==17|dat.residence$f.20118.0.0==18,"Scotland", NA))
dat.followup <- merge(dat.followup , dat.residence[,c(1,3)], by.x="f.eid", by.y="f.eid")

dat.followup$endtime <- ifelse(dat.followup$redince == "England", as.Date("2020-02-29"),
                               ifelse(dat.followup$redince =="Scotland", as.Date("2021-01-31"), 0))
dat.followup$endtime <- as.Date(dat.followup$endtime, "%Y-%m-%d",origin="1970-01-01")
dat.followup$endtime[is.na(dat.followup$endtime)] <- as.Date("2021-04-05") # the latest date of the CVD followup date
dat.followup$endtime_difftime <- as.numeric(difftime(dat.followup$endtime , dat.followup$entertime,units="days"), units="days")


write.csv(dat.followup, "C:/Users/sunx3/OneDrive - Memorial Sloan Kettering Cancer Center/My Work/Project/Else/RA and CVD_genetic correlation/data/UKBB/ukbb_RA&CAD_followup_2.csv",row.names = F)



