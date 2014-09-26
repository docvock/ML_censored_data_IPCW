rpart_functions <- function(varnms, cens_method = "IPCW", train.dat, test.dat = NULL,...) {  

	fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
	if(cens_method == "IPCW") {
		RP <- rpart(fmla, data = train.dat, weights = wts, method="class", control=rpart.control(...))
	}
	if(cens_method == "ZERO") {
		RP <- rpart(fmla, data = train.dat, method = "class", control = rpart.control(...))
	}
	if(cens_method == "DISCARD") {
		RP <- rpart(fmla, data = filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0)), method = "class",
			control=rpart.control(...))
	}	
  if(is.null(test.dat)) {
    probs <- predict(RP, type = "prob", na.action = na.omit)[,1]    
  } else {
    probs <- predict(RP, newdata = test.dat, type = "prob", na.action = na.omit)[,1]
  }
  return(probs)
}

## IPCW logisitc regression
logistic_functions <- function(varnms, cens_method = "IPCW", train.dat, test.dat = NULL, ID = 1:nrow(train.dat),...) {  

  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
  
	if (cens_method == "IPCW") {
		mydesign <- svydesign(id = ~ID, weights = ~wts, data = train.dat) 
		LR <- svyglm(fmla, family = quasibinomial(), design = mydesign) 
	}
	if (cens_method == "ZERO") {
		LR <- glm(fmla,data=train.dat,family="binomial")
	}
	if (cens_method == "DISCARD") {
		LR <- glm(fmla, data = filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0)),family = "binomial")
	}
	
  if(is.null(test.dat)) {
    probs <- predict(LR, type = "response", na.action=na.omit)    
  } else {
    probs <- predict(LR, newdata=test.dat, type = "response", na.action = na.omit)    
  }
  return(1-probs)
}

## Define the calibration statistic
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

Cindex <- function(p, test.dat) {
	test.dat$p <- p
	Cindex_all <- survConcordance(Surv(T.use, C.use) ~ p, data = test.dat)
	Cindex <- Cindex_all$concordance
	return(Cindex)
}
	

