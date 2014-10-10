code.base <- "C:/Users/Julian/Google Drive/iPredict/manuscripts/ML-for-censored-data/code/"
dat.base <- "H:/Carlson/Portal/Data/DataLock/DataLockv2.3.1/"

library(data.table)
## Read the data
HP.orig <- data.frame(fread(paste0(dat.base,"Data_ver-2.3.1_updated.csv"),verbose=TRUE))
HP.orig.sub <- subset(HP.orig,X.age>=40 & X.Comorbitity_All==0) ## Only those: Over age 40, no comorbidities

summary(HP.orig.sub)

HP.imp <- fread(paste0(dat.base,"Data_ver-2.3.1_imputed.csv"),verbose=TRUE)
HP <- subset(HP.imp,age>=40 & Comorbidity_All==0) ## Only those: Over age 40, no comorbidities


T.all <- HP$DaysToEvent_Fram
C.all <- HP$CVDEvent_Fram

postscript(paste0(code.base,"../FollowUp.eps"),width=6,height=6,paper="special",horizontal=FALSE,onefile=FALSE)
### Create a histogram of follow-up times
hist(T.all/365,xlab="Years of Follow-Up",ylab="Number of Patients",main=NULL,col="gray",breaks=10)
dev.off()
### Plot the percentage of patients with unknown event status
postscript(paste0(code.base,"../UnknownEvents.eps"),width=6,height=6,paper="special",horizontal=FALSE,onefile=FALSE)
times <- seq(0,max(T.all),length.out=100)
pct.unknown <- sapply(times,function(t) { mean(T.all<t & C.all==0)})
plot(times/365,pct.unknown,type="l",ylab="Percentage of patients with unknown event status",xlab="Time horizon (in years)")
dev.off()

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
