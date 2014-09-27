############################
#### Fit the HP dataset ####
############################

setwd("//hifs00.ahc.umn.edu/Data/Carlson/Portal")
dat.base <- "Data/DataLock/DataLockv2.3.1/"
source("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\code\\IPCW_fns_complete.R")

# libraries used
library(data.table)
library(dplyr) 
library(survey) # for weighted generalized linear models
library(rpart) # for recursive partitioning
library(earth) # for multivariate addaptive regression splines
library(mgcv) # for generalized additive models
library(caret) # for k-nearest neighbors
library(survival) # for Kaplan-Meier estimators

#  Read the data
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
#  Restrict Dataset to those over 40 and with no comorbidites
HP <- subset(HP.imp, age >= 40 & Comorbidity_All == 0) 
#  Define factor variables, get standardized variabels, remove outlier for Total Cholesterol 
HP <- mutate(HP, TC = pmin(TC, 500), gender_f = as.factor(gender),
	Has_Diab_f = as.factor(Has_Diab),
	Smoking_f = as.factor(Smoking),
	SBP_Meds_f = as.factor(SBP_Meds),
	gender_SBP_Meds = ifelse(gender==1 & SBP_Meds == 1, 1,
															ifelse(gender == 0 & SBP_Meds == 1, 2,
																ifelse(gender == 1 & SBP_Meds == 0, 3, 4))),
	gender_s = stand(gender), 
	age_s = stand(age),
	SBP_s = stand(SBP),
	BMI_s = stand(BMI),
	HDL_s = stand(HDL),
	TC_s = stand(TC),
	Has_Diab_s = stand(Has_Diab),
	Smoking_s = stand(Smoking),
	SBP_Meds_s = stand(SBP_Meds))

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
varnms <- c("gender_f", "age", "SBP", "BMI", "HDL", "TC", "Has_Diab_f", "Smoking_f", "SBP_Meds_f")
varnms_a <- c("gender_f", "age", "SBP", "BMI", "HDL", "TC", "Has_Diab_f", "Smoking", "SBP_Meds_f")
intx.varnms <- c(varnms,paste0("gender_f:",varnms[-1]),"SBP:SBP_Meds_f","gender_f:SBP:SBP_Meds_f")
varnms_fixed <- c("gender_f", "Has_Diab_f", "Smoking_f", "SBP_Meds_f", "gender_f:Has_Diab_f", 
	"gender_f:Smoking_f", "gender_f:SBP_Meds_f")
varnms_smooth <- cbind(c("age", "SBP", "BMI", "HDL", "TC"),
	c("gender_f", "gender_SBP_Meds", "gender_f", "gender_f", "gender_f"))
varnms_stand <- c("gender_s", "age_s", "SBP_s", "BMI_s", "HDL_s", "TC_s", "Has_Diab_s", 
	"Smoking_s", "SBP_Meds_s")

# Get IPC Weights
sf.C <- survfit(Surv(T.use, 1-C.use) ~ 1, data = HP.train)
KM.C <- function(t) {
    KM.C <- sapply(1:length(t), function(pos) {
    	sf.C$surv[min(which(sf.C$time > t[pos]))]
    })
  }
HP.train <- mutate(HP.train, wts = ifelse(T.use < (SURVTIME*365) & C.use==0, 0,
	                           1 / KM.C(pmin(T.use , SURVTIME*365)) ))


#  Recursive Partitionin Models
lossmat <- rbind(c(0,1),c(5,0))
pRP <- rpart_functions(varnms_a, "IPCW", HP.train, HP.test, lossmat, "information", cp = 0.00001, 
	minbucket = 100)
pRP.zero <- rpart_functions(varnms_a, "ZERO", HP.train, HP.test, lossmat, "information", 
	cp = 0.00001, minbucket = 100)
pRP.disc <- rpart_functions(varnms_a, "DISCARD", HP.train, HP.test, lossmat, "information", 
	cp = 0.00001, minbucket = 100)

#  Bagged Recursive Partitionin Models
lossmat <- rbind(c(0,1),c(5,0))
set.seed(1101985)
pRP.bag <- rpart_bag_functions(varnms_a, 50, "IPCW", HP.train, HP.test, lossmat, "information",
	cp = 0.00001, minbucket = 100)
pRP.bag.zero <- rpart_bag_functions(varnms_a, 50, "ZERO", HP.train, HP.test, lossmat, "information", 
	cp = 0.00001, minbucket = 100)
pRP.bag.disc <- rpart_functions(varnms_a, 50, "DISCARD", HP.train, HP.test, lossmat, "information", 
	cp = 0.00001, minbucket = 100)

#  Logistic Regression Models
pLR <- logistic_functions(intx.varnms, "IPCW", HP.train, HP.test)
pLR.zero <- logistic_functions(intx.varnms, "ZERO", HP.train, HP.test)
pLR.disc <- logistic_functions(intx.varnms, "DISCARD", HP.train, HP.test)

#  Bagged Logistic Regression Models
set.seed(1101985)
pLR.bag <- logistic_bag_functions(intx.varnms, 50, "IPCW", HP.train, HP.test)
pLR.bag.zero <- logistic_bag_functions(intx.varnms, 50, "ZERO", HP.train, HP.test)
pLR.bag.disc <- logistic_bag_functions(intx.varnms, 50, "DISCARD", HP.train, HP.test)

#  Generalized Additive Models
pGAM <- GAM_functions(varnms_fixed, varnms_smooth, "IPCW", HP.train, HP.test)
pGAM.zero <- GAM_functions(intx.varnms, "ZERO", HP.train, HP.test)
pGAM.disc <- GAM_functions(intx.varnms, "DISCARD", HP.train, HP.test)

#  Bagged Generalized Additive Models
pGAM.bag <- GAM_bag_functions(varnms_fixed, varnms_smooth, 10, "IPCW", HP.train, 
	HP.test)

#  Bagged MARS Models
set.seed(1101985)
pMARS <- MARS_bag_functions(varnms, 50, "IPCW", HP.train, HP.test)
pMARS.zero <- MARS_bag_functions(varnms, 50, "ZERO", HP.train, HP.test)
pMARS.disc <- MARS_bag_functions(varnms, 50, "DISC", HP.train, HP.test)

#  Bagged knn
set.seed(1101985)
pknn <- knn_bag_functions(varnms_stand, 50, "IPCW", HP.train, HP.test, k = 100)
pknn.zero <- knn_bag_functions(varnms_stand, 50, "ZERO", HP.train, HP.test, k = 100)
pknn.disc <- knn_bag_functions(varnms_stand, 50, "DISC", HP.train, HP.test, k = 100)


#  Calibration statistics
risk.cutpts = c(Inf,1-c(0.05,0.1,0.15,0.2),-Inf)
library(survMisc)

calib.stat(pRP, risk.cutpts, HP.test)
calib.stat(pRP.zero, risk.cutpts, HP.test)
calib.stat(pRP.disc, risk.cutpts, HP.test)

calib.stat(pRP.bag, risk.cutpts, HP.test)
calib.stat(pRP.bag.zero, risk.cutpts, HP.test)
calib.stat(pRP.bag.disc, risk.cutpts, HP.test)

calib.stat(pLR, risk.cutpts, HP.test)
calib.stat(pLR.zero, risk.cutpts, HP.test)
calib.stat(pLR.disc, risk.cutpts, HP.test)

calib.stat(pLR.bag, risk.cutpts, HP.test)
calib.stat(pLR.bag.zero, risk.cutpts, HP.test)
calib.stat(pLR.bag.disc, risk.cutpts, HP.test)

calib.stat(pGAM, risk.cutpts, HP.test)
calib.stat(pGAM.zero, risk.cutpts, HP.test)
calib.stat(pGAM.discard, risk.cutpts, HP.test)

calib.stat(pGAM.bag, risk.cutpts, HP.test)
calib.stat(pMARS, risk.cutpts, HP.test)
calib.stat(pknn, risk.cutpts, HP.test)

# Evaluate the C-index
Cindex(pRP, HP.test)
Cindex(pRP.zero, HP.test)
Cindex(pRP.disc, HP.test)

Cindex(pRP.bag, HP.test)
Cindex(pRP.bag.zero, HP.test)
Cindex(pRP.bag.disc, HP.test)

Cindex(pLR, HP.test)
Cindex(pLR.zero, HP.test)
Cindex(pLR.disc, HP.test)

Cindex(pLR.bag, HP.test)
Cindex(pLR.bag.zero, HP.test)
Cindex(pLR.bag.disc, HP.test)

Cindex(pGAM, HP.test)
Cindex(pGAM.zero, HP.test)
Cindex(pGAM.disc, HP.test)

Cindex(pGAM.bag, HP.test)
Cindex(pMARS, HP.test)
Cindex(pknn, HP.test)

## Evaluate the reclassification performance of well-calibrated models


