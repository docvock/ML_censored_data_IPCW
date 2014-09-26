############################
#### Fit the HP dataset ####
############################

setwd("//hifs00.ahc.umn.edu/Data/Carlson/Portal")
dat.base <- "Data/DataLock/DataLockv2.3.1/"
source("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\code\\IPCW_fns_complete.R")

# libraries used
library(data.table)
library(dplyr)
#library(peperr)
library(survey) 
library(rpart)
library(gam)
library(earth)
library(survival)

#  Read the data
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
#  Restrict Dataset to those over 40 and with no comorbidites
HP <- subset(HP.imp, age >= 40 & Comorbidity_All == 0) 
#  Define factor variables and remove outlier for Total Cholesterol 
HP <- mutate(HP, TC = pmin(TC, 500), gender_f = as.factor(gender),
	Has_Diab_f = as.factor(Has_Diab),
	Smoking_f = as.factor(Smoking),
	SBP_Meds_f = as.factor(SBP_Meds),
	gender_SBP_Meds = ifelse(gender==1 & SBP_Meds == 1, 1,
															ifelse(gender == 0 & SBP_Meds == 1, 2,
																ifelse(gender == 1 & SBP_Meds == 0, 3, 4))),
																gender_SBP_Meds = as.factor(gender_SBP_Meds))
# Time horizon for prediction in years
SURVTIME <- 7 
# Define the outcome of interest (time to Hard CVD, Soft CVD, etc.) 
HP <- mutate(HP, T.use = DaysToEvent_Fram, C.use = CVDEvent_Fram, 
	eventsurv.ind = as.factor(C.use == 1 & T.use <= (SURVTIME*365)))

#set.seed(1101985)
set.seed(1281985)
frac.train <- 0.75
train.set <- sample(1:nrow(HP), floor(frac.train*nrow(HP)), replace=FALSE)
test.set <- setdiff(1:nrow(HP), train.set)

HP.train <- data.frame(HP[train.set,])
HP.test <- data.frame(HP[test.set,])

#  Features to use in analysis
varnms <- c("gender_f","age","SBP","BMI","HDL","TC","Has_Diab_f","Smoking_f", "SBP_Meds_f")
intx.varnms <- c(varnms,paste0("gender_f:",varnms[-1]))

# Get IPC Weights
sf.C <- survfit(Surv(T.use, 1-C.use) ~ 1, data = HP.train)
KM.C <- function(t) {
    sf.C$surv[min(which(sf.C$time > t))]
  }
HP.train <- mutate(HP.train, wts = ifelse(T.use < (SURVTIME*365) & C.use==0, 0,
	                           1 / KM.C(min(T.use , SURVTIME*365)) ))


#  Recursive Partitionin Models
pRP <- rpart_functions(varnms, "IPCW", HP.train, HP.test, cp = 0.00001, minbucket = 100)
pRP.zero <- rpart_functions(varnms, "ZERO", HP.train, HP.test, cp = 0.00001, minbucket = 100)
pRP.disc <- rpart_functions(varnms, "DISCARD", HP.train, HP.test, cp = 0.00001, minbucket = 100)

#  Logistic Regression Models
pLR <- logistic_functions(intx.varnms, HP.train, HP.test)
pLR.zero <- logistic_functions(intx.varnms,HP.train, HP.test)
pLR.disc <- logistic_functions(intx.varnms,HP.train, HP.test)

#  Generalized Additive Models

#  Bagged MARS Models

#  Bagged Recursive Partition Models


risk.cutpts = c(Inf,1-c(0.05,0.1,0.15,0.2),-Inf)
library(survMisc)

calib.stat(pRP, risk.cutpts, HP.test)
calib.stat(pRP.zero,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pRP.disc,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR.zero,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR.disc,risk.cutpts,SURVTIME*365,T.test,C.test)

## Evaluate the reclassification performance of well-calibrated models
