library(rpart)

IPCW_rpart <- function(T,C,varnms,tau,train.dat,test.dat=NULL,parmlist=NULL,...) {  
  sf <- survfit(Surv(T,C)~1)  
  sf.C <- survfit(Surv(T,1-C)~1)
  
  KM <- function(t) {
    sf$surv[min(which(sf$time > t))]
  }
  
  KM.C <- function(t) {
    sf.C$surv[min(which(sf.C$time > t))]
  }
  
  wts <- sapply(1:length(T),function(j) {
    if((T[j] < (tau) & C[j]==0) ) { return(0) }
    else { return(1/KM.C(min(T[j],tau)))}
  })
  
  excl <- (C==0 & T <= tau)
  eventsurv.ind <- as.factor(C==1 & T <= tau)
  train.dat$eventsurv <- eventsurv.ind
  fmla <- as.formula(paste0("eventsurv~",paste0(varnms,collapse="+")))
  lossmat <- rbind(c(0,1),c(10,0))
  RP <- rpart(fmla,data=subset(train.dat,!excl),weights = wts[!excl],method="class",parms=parmlist,control=rpart.control(...))
  
  if(is.null(test.dat)) {
    probs <- predict(RP,type="prob",na.action=na.omit)[,1]    
  } else {
    probs <- predict(RP,newdata=test.dat,type="prob",na.action=na.omit)[,1]
  }
  return(probs)
}

## Decision tree setting observations censored prior to tau to zero
ZERO_rpart <- function(T,C,varnms,tau,train.dat,test.dat=NULL,parmlist=NULL,...) {  

  eventsurv.ind <- as.factor(C==1 & T <= tau)
  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))

  lossmat <- rbind(c(0,1),c(10,0))
  RP <- rpart(fmla,data=train.dat,method="class",parms=parmlist,control=rpart.control(...))
  
  if(is.null(test.dat)) {
    probs <- predict(RP,type="prob",na.action=na.omit)[,1]    
  } else {
    probs <- predict(RP,newdata=test.dat,type="prob",na.action=na.omit)[,1]
  }
  return(probs)
}

## Decision tree discarding all "unknown" results
DISCARD_rpart <- function(T,C,varnms,tau,train.dat,test.dat=NULL,parmlist=NULL,...) {  
  
  excl <- (C==0 & T <= tau)
  eventsurv.ind <- as.factor(C==1 & T <= tau)
  train.dat$eventsurv <- eventsurv.ind
  fmla <- as.formula(paste0("eventsurv~",paste0(varnms,collapse="+")))
  lossmat <- rbind(c(0,1),c(10,0))
  RP <- rpart(fmla,data=subset(train.dat,!excl),method="class",parms=parmlist,control=rpart.control(...))
    
  if(is.null(test.dat)) {
    probs <- predict(RP,type="prob",na.action=na.omit)[,1]    
  } else {
    probs <- predict(RP,newdata=test.dat,type="prob",na.action=na.omit)[,1]
  }
  return(probs)
}

library(survey) 

## IPCW logistic regression
IPCW_logistic <- function(T,C,varnms,tau,train.dat,test.dat=NULL,ID=1:nrow(train.dat),...) {  

  sf <- survfit(Surv(T,C)~1)  
  sf.C <- survfit(Surv(T,1-C)~1)
  
  KM <- function(t) {
    sf$surv[min(which(sf$time > t))]
  }
  
  KM.C <- function(t) {
    sf.C$surv[min(which(sf.C$time > t))]
  }
  
  wts <- sapply(1:length(T),function(j) {
    if((T[j] < (tau) & C[j]==0) ) { return(0) }
    else { return(1/KM.C(min(T[j],tau)))}
  })
  
  eventsurv.ind <- as.factor(C==1 & T <= tau)
  
  train.dat$eventsurv <- eventsurv.ind
  train.dat$ID <- ID
  
  fmla <- as.formula(paste0("eventsurv~",paste0(varnms,collapse="+")))
    
  mydesign <- svydesign(id=~ID, weights=~wts, data=train.dat) 
  LR <- svyglm( fmla, family=quasibinomial(), design=mydesign) 
    
  if(is.null(test.dat)) {
    probs <- predict(LR,type="response",na.action=na.omit)    
  } else {
    probs <- predict(LR,newdata=test.dat,type="response",na.action=na.omit)    
  }
  return(1-probs)
}

## Logistic regression setting censored prior to tau to zero
ZERO_logistic <- function(T,C,varnms,tau,train.dat,test.dat=NULL,...) {  

  eventsurv.ind <- as.integer(C==1 & T <= tau)
  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
  
  LR <- glm(fmla,data=train.dat,family="binomial")
  
  if(is.null(test.dat)) {
    probs <- predict(LR,type="response",na.action=na.omit)    
  } else {
    probs <- predict(LR,newdata=test.dat,type="response",na.action=na.omit)    
  }
  return(1-probs)
}

## Logistic regression discarding censored prior to tau
DISCARD_logistic <- function(T,C,varnms,tau,train.dat,test.dat=NULL,...) {  
  
  eventsurv.ind <- as.factor(C==1 & T <= tau)
  
  train.dat$eventsurv <- eventsurv.ind
  fmla <- as.formula(paste0("eventsurv~",paste0(varnms,collapse="+")))
  
  LR <- glm(fmla,data=subset(train.dat,!(T<=tau & C==0)),family="binomial")
  
  if(is.null(test.dat)) {
    probs <- predict(LR,type="response",na.action=na.omit)    
  } else {
    probs <- predict(LR,newdata=test.dat,type="response",na.action=na.omit)    
  }
  return(1-probs)
}

COX <- function(T,C,varnms,tau,train.dat,test.dat=NULL,...) {
  fmla <- as.formula(paste0("Surv(T,C)~",paste0(varnms,collapse="+")))
  Cox.train <- coxph(fmla,data=train.dat)
  p1 <- survfit(Cox.train,test.dat[1,varnms])
  ind.surv <- max(which(p1$time<=SURVTIME*365))

  if(is.null(test.dat)) {
    probs <- survfit(Cox.train,se.fit=FALSE,conf.type="none")$surv[ind.surv,]  
  } else {
    probs <- survfit(Cox.train,test.dat[,varnms],se.fit=FALSE,conf.type="none")$surv[ind.surv,]
  }
  return(probs)
}
