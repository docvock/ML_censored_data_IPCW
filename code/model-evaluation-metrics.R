### Function for computing risk reclassification - version which doesn't account for censoring
NRI <- function(risks.1,risks.2,event.ind) {
  event.reclass <- (sum(risks.2[event.ind]>risks.1[event.ind]) - sum(risks.1[event.ind] > risks.2[event.ind]))/sum(event.ind)
  nevent.reclass <- (sum(risks.2[!event.ind]<risks.1[!event.ind]) - sum(risks.1[!event.ind] < risks.2[!event.ind]))/sum(!event.ind)
  c(NRI.events=event.reclass,NRI.nonevents=nevent.reclass,NRI=event.reclass + nevent.reclass)
}

## Function for computing net reclassification improvement, taking censoring into account
## Due to Pencina (2011)
cNRI <- function(p1,p2,cutpts,test.dat,t) {
  risks.1 <- cut(1-p1,breaks=cutpts,labels=FALSE)
  risks.2 <- cut(1-p2,breaks=cutpts,labels=FALSE)
  T <- test.dat$T.use
  C <- test.dat$C.use
  
  n <- length(T)

  up.class <- (risks.2>risks.1)
  down.class <- (risks.2<risks.1)
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

## Harrell's C-index
Cindex <- function(p, test.dat) {
  test.dat$p <- 1-p
  Cindex_all <- survConcordance(Surv(T.use, C.use) ~ p, data = test.dat)
  Cindex <- Cindex_all$concordance
  return(Cindex)
}

## Calibration statistic
calib.stat <- function(p, cutpts, test.dat) {
  risk.class <- cut(p, cutpts, labels=FALSE)
  lev.stats <- sapply(1 : (length(cutpts) - 1), function(f) {
    ind <- which(risk.class == f)
    S.KM <- calcSurv(Surv(test.dat$T.use[ind], test.dat$C.use[ind]))
    ind.surv <- max(which(S.KM$t <= (SURVTIME*365)))
    p.KM <- S.KM$SKM[ind.surv]
    ##print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv])))
    (mean(p[ind], na.rm = TRUE) - p.KM)^2 / S.KM$SKMV[ind.surv]
  })
  
  calib.stat <- (sum(lev.stats, na.rm=TRUE))
  return(calib.stat)
}

obsPred <- function(p,cutpts,t) {
  risk.class <- cut(p,cutpts,labels=FALSE)
  DF <- data.frame(t(sapply(1:(length(cutpts)-1),function(f) {
    ind <- which(risk.class==f)
    S.KM <- calcSurv(Surv(T.test[ind],C.test[ind]))
    ind.surv <- max(which(S.KM$t<=t))
    p.KM <- S.KM$SKM[ind.surv]
    mean.pred <- mean(p[ind],na.rm=TRUE)
    c(p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv]))
  })))
  colnames(DF) <- c("obs","pred","var")
  DF
}