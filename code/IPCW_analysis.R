############################
#### Fit the HP dataset ####
############################

setwd("//hifs00.ahc.umn.edu/Data/Carlson/Portal")
dat.base <- "Data/DataLock/DataLockv2.3.1/"
source("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\IPCW_fns.R")
#source("C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/IPCW_fns_complete.R")
#source("C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/model-evaluation-metrics.R")
source("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\code\\model-evaluation-metrics.R")

# libraries used
library(data.table)
library(dplyr) 
library(survey) # for weighted generalized linear models
library(rpart) # for recursive partitioning
library(mgcv) # for generalized additive models
library(yaImpute) # for knn estimators
library(survival) # for Kaplan-Meier estimators
library(xtable) # create tables for latex
library(caret)
library(mvtnorm) 
library(survMisc)
quantile <- stats::quantile
cancor <- stats::cancor

#  Read the data
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
#  Restrict Dataset to those over 40 and with no comorbidites
HP <- subset(HP.imp, age >= 40 & Comorbidity_All == 0) 
#  Define factor variables, get standardized variabels, remove outlier for Total Cholesterol 
HP <- mutate(HP, 
	#TC = pmin(TC, 500), 
	gender_f = as.factor(cut(gender, breaks = c(-1,0,1), labels = FALSE)),
	Has_Diab_f = as.factor(cut(Has_Diab, breaks = c(-1,0,1), labels = FALSE)),
	Smoking_f = as.factor(cut(Smoking, breaks = c(-1,0,1,2), labels = FALSE)),
	SBP_Meds_f = as.factor(cut(SBP_Meds, breaks = c(-1,0,1), labels = FALSE)),
	BMI_f = cut(BMI, breaks = c(0, 25, 30, 35, 40, 100), labels = FALSE),
	age_f = cut(age, breaks = c(0, 50, 60, 70, 80, 110), labels = FALSE),
	SBP = pmin(SBP, quantile(SBP, p=c(0.995))),
	HDL = pmin(HDL, quantile(HDL, p=c(0.995))),
	TC = pmin(TC, quantile(TC, p=c(0.995))),
	log_SBP = log(SBP),
	log_HDL = log(HDL),
	log_TC = log(TC),	
	Smoking_bin = ifelse(Smoking == 1, 1, 0)
	)

# Time horizon for prediction in years
SURVTIME <- 5 
# Define the outcome of interest (time to Hard CVD, Soft CVD, etc.) 
HP <- mutate(HP, T.use = DaysToEvent_Fram, C.use = CVDEvent_Fram, 
	eventsurv.ind = as.factor(C.use == 1 & T.use <= (SURVTIME*365)))

set.seed(1281985)
frac.train <- 0.75
train.set <- sample(1:nrow(HP), floor(frac.train*nrow(HP)), replace=FALSE)
test.set <- setdiff(1:nrow(HP), train.set)

HP.train <- data.frame(HP[train.set,])
HP.test <- data.frame(HP[test.set,])

#  Features to use in analysis
varnms <- c("gender_f", "age", "SBP", "BMI", "HDL", "TC", "Has_Diab_f", "Smoking_f", "SBP_Meds_f")
varnms_fixed <- c("gender_f", "Has_Diab_f", "Smoking_f", "SBP_Meds_f")
varnms_smooth <- c("age", "SBP", "BMI", "HDL", "TC")
varnms_stand <- c("gender", "age", "SBP", "BMI", "HDL", "TC", "SBP_Meds", "Has_Diab",
	"Smoking_bin")

# Get IPC Weights
sf.C <- survfit(Surv(T.use, 1-C.use) ~ 1, data = HP.train)
KM.C <- function(t) {
    KM.C <- sapply(1:length(t), function(pos) {
    	sf.C$surv[min(which(sf.C$time > t[pos]))]
    })
  }
HP.train <- mutate(HP.train, wts = ifelse(T.use < (SURVTIME*365) & C.use==0, 0,
	                           1 / KM.C(pmin(T.use , SURVTIME*365)) ))

# Risk Cut-Offs (used in CV)
risk.cutpts = c(Inf,1-c(0.05, 0.1, 0.15, 0.20),-Inf)

#  Time needed to train all models
start <- proc.time()

#  Recursive Partitioning Models
pRP <- rpart_functions(varnms, "IPCW", train.dat = HP.train, test.dat = HP.test, myseed = 1101985,
                       losscost_seq = c(2.5, 5, 10), minbucket_seq = c(50, 100, 200), 
                       split = "information", risk.cutpts = risk.cutpts, cp=1e-4)
pRP.zero <- rpart_functions(varnms, "ZERO", train.dat = HP.train, test.dat = HP.test, myseed = 1101985,
                       losscost_seq = c(2.5, 5, 10), minbucket_seq = c(50, 100, 200), 
                       split = "information", risk.cutpts = risk.cutpts, cp=1e-4)
pRP.disc <- rpart_functions(varnms, "DISCARD", train.dat = HP.train, test.dat = HP.test, myseed = 1101985,
                       losscost_seq = c(2.5, 5, 10), minbucket_seq = c(50, 100, 200), 
                       split = "information", risk.cutpts = risk.cutpts, cp=1e-4)
pRP.split <- rpart_functions(varnms, "SPLIT", train.dat = HP.train, test.dat = HP.test, myseed = 1101985,
                       losscost_seq = c(2.5, 5, 10), minbucket_seq = c(50, 100, 200), 
                       split = "information", risk.cutpts = risk.cutpts, cp=1e-4)

#  Logistic Regression Models
pLR <- as.vector(logistic_functions(varnms, "IPCW", HP.train, HP.test))
pLR.zero <- as.vector(logistic_functions(varnms, "ZERO", HP.train, HP.test))
pLR.disc <- as.vector(logistic_functions(varnms, "DISCARD", HP.train, HP.test))
pLR.split <- as.vector(logistic_functions(varnms, "SPLIT", HP.train, HP.test))

#  Generalized Additive Models
pGAM <- GAM_functions(varnms_fixed, varnms_smooth, "IPCW", HP.train, HP.test)
pGAM.zero <- GAM_functions(varnms_fixed, varnms_smooth, "ZERO", HP.train, HP.test)
pGAM.disc <- GAM_functions(varnms_fixed, varnms_smooth, "DISCARD", HP.train, HP.test)
pGAM.split <- GAM_functions(varnms_fixed, varnms_smooth, "SPLIT", HP.train, HP.test)

#  kNN Models
pkNN <- knn_functions(varnms_stand, kmax = 1000, cens_method = "IPCW", dist_method = "msn", 
                          train.dat = HP.train, test.dat = HP.test, risk.cutpts = risk.cutpts, 
													myseed = 1101985, n.folds = 5, tuning.rule = "1SE")
pkNN.zero <- knn_functions(varnms_stand, kmax = 1000, cens_method = "ZERO", dist_method = "msn", 
                          train.dat = HP.train, test.dat = HP.test, risk.cutpts = risk.cutpts, 
													myseed = 1101985, n.folds = 5, tuning.rule = "1SE")
pkNN.disc <- knn_functions(varnms_stand, kmax = 1000, cens_method = "DISCARD", dist_method = "msn", 
                          train.dat = HP.train, test.dat = HP.test, risk.cutpts = risk.cutpts, 
													myseed = 1101985, n.folds = 5, tuning.rule = "1SE")
pkNN.split <- knn_functions(varnms_stand, kmax = 1000, cens_method = "SPLIT", dist_method = "msn", 
                          train.dat = HP.train, test.dat = HP.test, risk.cutpts = risk.cutpts, 
													myseed = 1101985, n.folds = 5, tuning.rule = "1SE")

#  Bayesian Network Models
pBayes <- Bayes_functions("IPCW", HP.train, HP.test)
pBayes.zero <- Bayes_functions("ZERO", HP.train, HP.test)
pBayes.disc <- Bayes_functions("DISCARD", HP.train, HP.test)
pBayes.split <- Bayes_functions("SPLIT", HP.train, HP.test)

# Cox model
fmla <- as.formula(paste0("Surv(T.use,C.use)~",paste0(varnms,collapse="+")))
Cox.train <- coxph(fmla,data=HP.train)
p1 <- survfit(Cox.train,HP.test[1,varnms])
ind.surv <- max(which(p1$time<=SURVTIME*365))

pCPH <- survfit(Cox.train,HP.test[,varnms],se.fit=FALSE,conf.type="none")$surv[ind.surv,]

#  Time needed to train all models
proc.time()-start

#######  Compute calibration and reclassification statistics

#risk.cutpts = c(Inf,1-c(0.05, 0.1, 0.15, 0.20),-Inf)

ML.preds <- list(pRP, pRP.disc, pRP.zero, pRP.split, pkNN, pkNN.disc, pkNN.zero, pkNN.split, 
	pBayes, pBayes.disc, pBayes.zero, pBayes.split, 
	pLR, pLR.disc, pLR.zero, pLR.split, 
	pGAM, pGAM.disc, pGAM.zero, pGAM.split,
	pCPH)

pred.rate <- 100*(1-unlist(lapply(ML.preds,mean)))
calib <- unlist(lapply(ML.preds,calib.stat,cutpts=risk.cutpts,test.dat=HP.test))
c.index <- unlist(lapply(ML.preds,Cindex,test.dat=HP.test))

# Calibration and C-index Tables
results.tab <- data.frame(pred.rate, calib, c.index)
results.tab <- rbind(rep(NA, 3), results.tab[1:4, ], 
										rep(NA, 3), results.tab[5:8, ], 
										rep(NA, 3), results.tab[9:12, ], 
										rep(NA, 3), results.tab[13:16, ], 
										rep(NA, 3), results.tab[17:20, ], 
										rep(NA, 3), results.tab[21, ])
method <- c("\\hspace{0.5in} IPCW", "\\hspace{0.5in} Discard", "\\hspace{0.5in} Zero",  "\\hspace{0.5in} Split")
rownms <- c("\\bf{Tree}", method, "\\bf{k-NN", method, "\\bf{Bayes}", method, "\\bf{Logistic}", 
	method, "\\bf{GAM}", method, "\\bf{Cox}", "")

#rownms <- c(as.vector(outer(c("IPCW-", "Zero-", "Discard-", "Split-"), 
#	c("Tree", "Logistic", "GAM", "k-NN", "Bayes"), 
#	paste0)), "Cox")
results.tab <- data.frame(Method=rownms, results.tab)
colnms <- c("Method", "Predicted event rate (%)", "Calibration", "C-Index")
#rownames(results.tab) <- rownms
colnames(results.tab) <- colnms

print( xtable(results.tab,  
	align= "llccc",digits = c(0, 0, 2, 2, 3)),
  hline.after=c(0, 5, 10, 15, 20, 25), include.rownames = FALSE, sanitize.text.function = function(x){x})

## Calibration Plots

pdf("C:\\Users\\bstvock\\Documents\\research\\IPCW_machine_learning\\ML-for-censored-data\\calibration_plots.pdf",
	height = 15, width = 5)
par(mfrow=c(5,1))
# Recursive Partitioning
calib.plot(ML.preds[1:4], cutpts=risk.cutpts, test.dat=HP.test, 
	main.title = "Calibration of Classification Trees", pch.use = 16:19, 
	col.use = c("black", "orange", "green", "blue"), lty.use = 1:4,
	legend.use = TRUE, legend.text = c("IPCW", "Discard", "Zero", "Split"), legend.loc = "bottomleft")

# kNN Models
calib.plot(ML.preds[5:8], cutpts=risk.cutpts, test.dat=HP.test, 
	main.title = "Calibration of k-Neartest Neighbor", pch.use = 16:19, 
	col.use = c("black", "orange", "green", "blue"), lty.use = 1:4,
	legend.use = TRUE, legend.text = c("IPCW", "Discard", "Zero", "Split"))

# Bayesian Network Models
calib.plot(ML.preds[9:12], cutpts=risk.cutpts, test.dat=HP.test, 
	main.title = "Calibration of Bayesian Network Models", pch.use = 16:19, 
	col.use = c("black", "orange", "green", "blue"), lty.use = 1:4,
	legend.use = TRUE, legend.text = c("IPCW", "Discard", "Zero", "Split"))

# Logistic Regression
calib.plot(ML.preds[13:16], cutpts=risk.cutpts, test.dat=HP.test, 
	main.title = "Calibration of Logistric Regression Models", pch.use = 16:19, 
	col.use = c("black", "orange", "green", "blue"), lty.use = 1:4,
	legend.use = TRUE, legend.text = c("IPCW", "Discard", "Zero", "Split"))

# Generalized Additive Models
calib.plot(ML.preds[17:20], cutpts=risk.cutpts, test.dat=HP.test, 
	main.title = "Calibration of Generalized Additive Models", pch.use = 16:19, 
	col.use = c("black", "orange", "green", "blue"), lty.use = 1:4,
	legend.use = TRUE, legend.text = c("IPCW", "Discard", "Zero", "Split"))

dev.off()

### Compute the cNRI

predmat <- data.frame(Tree = pRP, kNN = pkNN, Bayes = pBayes, Logistic = pLR,  GAM = pGAM, Cox = pCPH)
c.nri <- -combn(1:ncol(predmat),2,function(x) { 
  cNRI(p1 = predmat[, x[[1]]], p2 = predmat[, x[[2]]], cutpts = 1-risk.cutpts, test.dat = HP.test,
  	t=SURVTIME*365, class_wt = c(1, 1))
})
c.nri.nms <- combn(1:ncol(predmat),2,function(x) { paste(colnames(predmat[,x]), collapse=" vs. ") })
colnames(c.nri) = c.nri.nms

### Compute the weighted cNRI

  KM <- survfit(Surv(T.use, C.use) ~ 1, data = HP.test) ## overall
  p.KM <- 1 - KM$surv[max(which(KM$time<=SURVTIME*365))] ## P(event)

weight.c.nri <- -combn(1:ncol(predmat),2,function(x) { 
  cNRI(p1 = predmat[, x[[1]]], p2 = predmat[, x[[2]]], cutpts = 1-risk.cutpts, test.dat = HP.test,
  	t=SURVTIME*365, class_wt = c(p.KM, 1 - p.KM))
})
c.nri.nms <- combn(1:ncol(predmat),2,function(x) { paste(colnames(predmat[,x]), collapse=" vs. ") })
colnames(weight.c.nri) = c.nri.nms

### Prepare results tables

### Net reclassification
colnms <- c("cNRI (Events)", "cNRI (Non-Events)", "cNRI (Overall)", "cNRI(Overall Weighted)")
results.tab <- cbind(t(c.nri), weight.c.nri[3, ])
colnames(results.tab) <- colnms

print(xtable(results.tab, align = "lcccc", digits = c(0,3,3,3,3)) )

### Weighted Net reclassification
results.tab <- t(weight.c.nri)
colnames(results.tab) <- colnms

print( xtable(results.tab, align = "lccc", digits = c(0, 3, 3, 3)) )
