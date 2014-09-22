############################
#### Fit the HP dataset ####
############################

code.base <- "C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/"
dat.base <- "H:/Carlson/Portal/Data/DataLock/DataLockv2.3.1/"
source(paste0(code.base,"IPCW-fns.R"))

library(data.table)
## Read the data
#HP <- fread(paste0(dat.base,"Data_ver-2.3.1_updated.csv"),verbose=TRUE)
HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
HP <- subset(HP.imp,age>=40&Comorbidity_All==0) ## Only those over age 40 without comorbidities

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

pRP <- IPCW_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,cp=0.00001,minbucket=100)
pRP.zero <- ZERO_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,cp=0.00001,minbucket=100)
pRP.disc <- DISCARD_rpart(T.train,C.train,varnms,SURVTIME*365,HP.train,HP.test,cp=0.00001,minbucket=100)

intx.varnms <- c(varnms,paste0("gender:",varnms[-1]))

pLR <- IPCW_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)
pLR.zero <- ZERO_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)
pLR.disc <- DISCARD_logistic(T.train,C.train,intx.varnms,SURVTIME*365,HP.train,HP.test)

risk.cutpts = c(Inf,1-c(0.05,0.1,0.15,0.2),-Inf)

calib.stat(pRP,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pRP.zero,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pRP.disc,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR.zero,risk.cutpts,SURVTIME*365,T.test,C.test)
calib.stat(pLR.disc,risk.cutpts,SURVTIME*365,T.test,C.test)

## Evaluate the reclassification performance of well-calibrated models
