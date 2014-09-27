############################
#### Fit the HP dataset ####
############################

code.base <- "C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/"
dat.base <- "H:/Carlson/Portal/Data/DataLock/DataLockv2.3.1/"
source(paste0(code.base,"calibration.R"))
source(paste0(code.base,"reclassification.R"))
source(paste0(code.base,"IPCW-fns.R"))

library(data.table)
## Read the data
#HP <- fread(paste0(dat.base,"Data_ver-2.3.1_updated.csv"),verbose=TRUE)
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
HP <- subset(HP.imp,age>=40 & Comorbidity_All==0) ## Only those: Over age 40, no comorbidities

HP$Smoking1 <- as.integer(HP$Smoking==1)
HP$Smoking2 <- as.integer(HP$Smoking==2)

varnms <- c("gender","age","SBP","BMI","HDL","TC","SBP_Meds","Smoking1","Smoking2","Has_Diab")

frac.train <- 0.75
train.set <- sample(1:nrow(HP),frac.train*nrow(HP),replace=FALSE)
test.set <- setdiff(1:nrow(HP),train.set)

HP.train <- data.frame(HP[train.set,])
HP.test <- data.frame(HP[test.set,])


T.train <- HP.train$DaysToEvent_Fram
C.train <- HP.train$CVDEvent_Fram

library(peperr)

T.test <- HP.test$DaysToEvent_Fram
C.test <- HP.test$CVDEvent_Fram

SURVTIME <- 5 ## Time horizon for prediction in years

### Fit the regression tree models
lossmat <- rbind(c(0,1),c(10,0))
pRP <- IPCW_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="gini"),
                  cp=0,minbucket=200)
pRP.zero <- ZERO_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="gini"),
                  cp=0,minbucket=200)
pRP.disc <- DISCARD_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="gini"),
                  cp=0,minbucket=200)

### Fit the logistic models
pLR <- IPCW_logistic(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test)
pLR.zero <- ZERO_logistic(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test)
pLR.disc <- DISCARD_logistic(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test)

### Fit a Cox model
pCOX <- COX(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test)

### Fit logistic regression with interaction terms
# intx.varnms <- c(varnms,paste0("gender:",varnms[-1]))
# 
# pLR <- IPCW_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)
# pLR.zero <- ZERO_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)
# pLR.disc <- DISCARD_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)

## Define a results table

results.tab <- data.frame(method=c("IPCW-Tree","Zero-Tree","Discard-Tree","IPCW-Logistic","Zero-Logistic","Discard-Logistic","Cox"))

### Compute the overall predicted event rate
event.rates <- 1-unlist(lapply(list(pRP,pRP.zero,pRP.disc,pLR,pLR.zero,pLR.disc,pCOX),mean))
results.tab$pct <- 100*event.rates

### Compute the calibration
risk.cutpts = c(Inf,1-c(0.05,0.1,0.15,0.2),-Inf)

calib.stats <- unlist(lapply(list(pRP,pRP.zero,pRP.disc,pLR,pLR.zero,pLR.disc,pCOX),function(prob){
  calib.stat(prob,risk.cutpts,SURVTIME*365,T.test,C.test,printTable=TRUE)
}))

results.tab$calib <- calib.stats

## Evaluate the NRI of well-calibrated models
brks = 1-risk.cutpts

rRP <- cut(1-pRP,breaks=brks,labels=FALSE)
rRP.zero <- cut(1-pRP.zero,breaks=brks,labels=FALSE)
rRP.disc <- cut(1-pRP.disc,breaks=brks,labels=FALSE)

rLR <- cut(1-pLR,breaks=brks,labels=FALSE)
rLR.zero <- cut(1-pLR.zero,breaks=brks,labels=FALSE)
rLR.disc <- cut(1-pLR.disc,breaks=brks,labels=FALSE)

rCOX <- cut(1-pCOX,breaks=brks,labels=FALSE)

cIndex <- lapply(list(pRP,pRP.zero,pRP.disc,pLR,pLR.zero,pLR.disc,pCOX),function(pred) {
  ind <- !is.na(pred)
  compute.cIndex(1-pred[ind],T.test[ind],C.test[ind],SURVTIME*365)
})

results.tab$cIndex <- unlist(cIndex)

cNRI <- lapply(list(rRP,rLR),function(rP) {
  ind <- !is.na(rP)
  compute.cNRI(rCOX[ind],rP[ind],T.test[ind],C.test[ind],SURVTIME*365)
})  

results.tab$cNRI <- rep(NA,nrow(results.tab))
results.tab$cNRI[c(1,4)] <- do.call(rbind,cNRI)[,3] 

print(results.tab)
summ <- summary(survfit(Surv(T.test,C.test)~1))
## Overall event probability
summ$surv[max(which(summ$time<=SURVTIME*365))]

print(xtable(results.tab,align="ll|cccc",digits=c(0,0,2,3,3,3)),include.rownames=FALSE,only.contents=TRUE)
