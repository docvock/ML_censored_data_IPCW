code.base <- "C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/"
dat.base <- "H:/Carlson/Portal/Data/DataLock/DataLockv2.3.1/"

library(data.table)
## Read the data
#HP <- fread(paste0(dat.base,"Data_ver-2.3.1_updated.csv"),verbose=TRUE)

HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
HP <- subset(HP.imp,age>=40 & Comorbidity_All==0) ## Only those: Over age 40, no comorbidities

frac.train <- 0.75
train.set <- sample(1:nrow(HP),frac.train*nrow(HP),replace=FALSE)
test.set <- setdiff(1:nrow(HP),train.set)

T.all <- HP$DaysToEvent_Fram
C.all <- HP$CVDEvent_Fram

HP.train <- data.frame(HP[train.set,])
HP.test <- data.frame(HP[test.set,])

T.train <- HP.train$DaysToEvent_Fram
C.train <- HP.train$CVDEvent_Fram

T.test <- HP.test$DaysToEvent_Fram
C.test <- HP.test$CVDEvent_Fram

## % censored prior to five years
mn.cens <- mean(C.train==0 & T.train <= (5*365))
print(mn.cens)

## Number evaluated in discard
sum.uncens <- sum(!(C.train==0 & T.train <= (5*365)))
print(sum.uncens)

### Get numbers to tabulate the predictors
summary(HP)
table(HP$gender)
table(HP$gender)/nrow(HP)
table(HP$Smoking)
table(HP$Smoking)/nrow(HP)
table(HP$SBP_Meds)
table(HP$SBP_Meds)/nrow(HP)
table(HP$Has_Diab)
table(HP$Has_Diab)/nrow(HP)

### Total number of events
sum(C.all==1&T.all<=(5*365))

## Overall KM event rate
sf <- survfit(Surv(T.all,C.all)~1)
print(1-sf$surv[max(which(sf$time<=(5*365)))])
