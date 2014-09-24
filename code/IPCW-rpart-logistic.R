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
HP <- subset(HP.imp,Comorbidity_All==0) ## Only those over age 40 without comorbidities

frac.train <- 0.75
train.set <- sample(1:nrow(HP),frac.train*nrow(HP),replace=FALSE)
test.set <- setdiff(1:nrow(HP),train.set)

HP.train <- data.frame(HP[train.set,])
HP.test <- data.frame(HP[test.set,])

## Fit the models
SURVTIME <- 7

T.train <- HP.train$DaysToEvent_Fram
C.train <- HP.train$CVDEvent_Fram

varnms <- c("gender","age","SBP","BMI","HDL","TC","Has_Diab","Smoking")
##varnms <- c("age","gender","SBP","BMI","HDL","TC") ## Reasonably well-calibrated

library(peperr)

T.test <- HP.test$DaysToEvent_Fram
C.test <- HP.test$CVDEvent_Fram

SURVTIME <- 7 ## Time horizon for prediction in years

### Fit the regression tree models
lossmat <- rbind(c(0,1),c(10,0))
pRP <- IPCW_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="information"),
                  cp=0.0001,minbucket=20)
pRP.zero <- ZERO_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="information"),
                  cp=0.0001,minbucket=20)
pRP.disc <- DISCARD_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,
                  parmlist=list(loss=lossmat,split="information"),
                  cp=0.0001,minbucket=20)

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

results.tab <- data.frame(method=c("rpart_ipcw","rpart_zero","rpart_discard","lr_ipcw","lr_zero","lr_disc","cox"))

### Compute the calibration

risk.cutpts = c(Inf,1-c(0.05,0.1,0.15,0.2),-Inf)

calib.stats <- unlist(lapply(list(pRP,pRP.zero,pRP.disc,pLR,pLR.zero,pLR.disc,pCOX),function(prob){
  calib.stat(prob,risk.cutpts,SURVTIME*365,T.test,C.test)
}))

results.tab$calib <- calib.stats

## Evaluate the reclassification performance of well-calibrated models
brks = 1-risk.cutpts

rRP <- cut(1-pRP,breaks=brks,labels=FALSE)
rLR <- cut(1-pLR,breaks=brks,labels=FALSE)
rCOX <- cut(1-pCOX,breaks=brks,labels=FALSE)

RR.c <- lapply(list(rRP,rLR),function(rP) {
  ind <- !is.na(rP)
  compute.cNRI(rCOX[ind],rP[ind],T.test[ind],C.test[ind],SURVTIME*365)
})  

results.tab$cNRI <- rep(NA,nrow(results.tab))
results.tab$cNRI[1] <- RR.c[[1]][3]
results.tab$cNRI[4] <- RR.c[[2]][3]


cIndex <- lapply(list(pRP,pLR,pCOX),function(pred) {
  ind <- !is.na(pred)
  compute.cIndex(pred[ind],T.test[ind],C.test[ind],SURVTIME*365)
})

results.tab$cIndex <- rep(NA,nrow(results.tab))
results.tab$cIndex[c(1,4,7)] <- unlist(cIndex)

print(results.tab)