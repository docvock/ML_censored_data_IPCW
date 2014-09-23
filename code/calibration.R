library(survMisc)
calib.stat <- function(p,cutpts,t) {
  risk.class <- cut(p,cutpts,labels=FALSE)
  lev.stats <- sapply(1:(length(cutpts)-1),function(f) {
    ind <- which(risk.class==f)
    S.KM <- calcSurv(Surv(T.test[ind],C.test[ind]))
    ind.surv <- max(which(S.KM$t<=t))
    p.KM <- S.KM$SKM[ind.surv]
    ##print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv])))
    (mean(p[ind],na.rm=TRUE) - p.KM)^2/S.KM$SKMV[ind.surv]
  })
  
  sum(lev.stats,na.rm=TRUE)
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