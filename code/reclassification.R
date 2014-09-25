### Function for computing risk reclassification
compute.NRI <- function(risks.1,risks.2,event.ind) {
  event.reclass <- (sum(risks.2[event.ind]>risks.1[event.ind]) - sum(risks.1[event.ind] > risks.2[event.ind]))/sum(event.ind)
  nevent.reclass <- (sum(risks.2[!event.ind]<risks.1[!event.ind]) - sum(risks.1[!event.ind] < risks.2[!event.ind]))/sum(!event.ind)
  c(NRI.events=event.reclass,NRI.nonevents=nevent.reclass,NRI=event.reclass + nevent.reclass)
}

## Function for computing net reclassification improvement, taking censoring into account
## Due to Wolfson
compute.wNRI <- function(risks.1,risks.2,T,C,t) {
  KM <- survfit(Surv(T,C)~1)
  p.fail.t <- 1 - KM$surv[max(which(KM$time <= t))]
  p.fail.obs <- sapply(T,function(t) {
    myind <- which(KM$time <= t)
    if(length(myind)==0) {
      out <- 0
      }
    else {
      out <- 1 - KM$surv[max(myind)]}
    out
  })
                       
  p.fail <- (p.fail.t - p.fail.obs) / (1 - p.fail.obs)
  ind <- which(C==0&T<t) ## individuals who are censored prior to t
  event.prob <- as.integer(C==1&T<t)
  event.prob[ind] <- p.fail[ind]
  
  nri.e <- ( sum((risks.2>risks.1)*event.prob) - sum((risks.1>risks.2)*event.prob) ) / sum(event.prob)
  nri.ne <- ( sum((risks.2<risks.1)*(1-event.prob)) - sum((risks.1<risks.2)*(1-event.prob)) ) / sum(1-event.prob)
  c(NRI.events=nri.e,NRI.nonevents=nri.ne,NRI=nri.e+nri.ne)
}

## Function for computing net reclassification improvement, taking censoring into account
## Due to Pencina (2011)
compute.cNRI <- function(risks.1,risks.2,T,C,t) {
  n <- length(T)

  up.class <- (risks.2>risks.1)
  down.class <- (risks.1>risks.2)
  n.u <- sum(up.class) ## Number up-classified
  n.d <- sum(down.class) ## Number down-classified

  KM <- survfit(Surv(T,C)~1) ## overall
  p.KM <- 1 - KM$surv[max(which(KM$time<=t))] ## P(event)
  
  if(n.u==0) {
    p.KM.u <- 1
    } else {
      KM.u <- survfit(Surv(T,C)~1,subset=up.class) ## up-classified
      p.KM.u <- 1 - KM.u$surv[max(which(KM.u$time<=t))] ## P(event|up)
    }
  if(n.d==0) {
    p.KM.d <- 1
  } else {
    KM.d <- survfit(Surv(T,C)~1,subset=down.class) ## down-classified
    p.KM.d <- 1 - KM.d$surv[max(which(KM.d$time<=t))] ## P(event|down)
  }

  nri.e <- (n.u*p.KM.u - n.d*p.KM.d)/(n*p.KM)
  nri.ne <- (n.d*(1-p.KM.d) - n.u*(1-p.KM.u))/(n*(1-p.KM))
  
  c(cNRI.events=nri.e,cNRI.nonevents=nri.ne,cNRI=nri.e+nri.ne)
}

## C-index for censored data due to Harrell
compute.cIndex <- function(pred,T,C,t) {
  
  ind1 <- (C==1&T<=t)
  ind2 <- ind1 | (T>t)
  A <- outer(T[ind1],T[ind2],FUN=function(x,y){as.integer(x<y)})
  B <- outer(pred[ind1],pred[ind2],FUN=function(x,y){as.integer(x<y)})
  
  num <- sum(A*B)
  denom <- length(A)
  num/denom
}
