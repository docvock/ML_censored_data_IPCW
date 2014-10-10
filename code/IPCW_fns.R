stand <- function(x){
	stand <- (x-mean(x))/sd(x)
}

rpart_functions <- function(varnms, cens_method = "IPCW", train.dat, test.dat = NULL, 
	lossmat = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), split = "information",...) {  

	fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
	if(cens_method == "IPCW") {
		RP <- rpart(fmla, data = train.dat, weights = wts, method="class", 
			parms = list(loss = lossmat, split = split), control=rpart.control(...))
	}
	if(cens_method == "ZERO") {
		RP <- rpart(fmla, data = train.dat, method = "class", 
			parms = list(loss = lossmat, split = split), control = rpart.control(...))
	}
	if(cens_method == "DISCARD") {
		RP <- rpart(fmla, data = filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0)), 
			method = "class", parms = list(loss = lossmat, split = split), control=rpart.control(...))
	}	
  if(is.null(test.dat)) {
    probs <- predict(RP, type = "prob", na.action = na.omit)[,1]    
  } else {
    probs <- predict(RP, newdata = test.dat, type = "prob", na.action = na.omit)[,1]
  }	
  return(probs)
}



## IPCW logistic regression
logistic_functions <- function(varnms, cens_method = "IPCW", train.dat, test.dat = NULL, 
	ID = 1:nrow(train.dat),...) {  

  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
  
	if (cens_method == "IPCW") {
		mydesign <- svydesign(id = ~ID, weights = ~wts, data = train.dat) 
		LR <- svyglm(fmla, family = quasibinomial(), design = mydesign) 
	}
	if (cens_method == "ZERO") {
		LR <- glm(fmla,data=train.dat,family="binomial")
	}
	if (cens_method == "DISCARD") {
		LR <- glm(fmla, data = filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0)), 
			family = "binomial")
	}
	
  if(is.null(test.dat)) {
    probs <- predict(LR, type = "response", na.action=na.omit)    
  } else {
    probs <- predict(LR, newdata=test.dat, type = "response", na.action = na.omit)
  }
  return(1-probs)
}


# IPCW Generalized Additive Model
GAM_functions <- function(varnms_fixed, varnms_smooth, cens_method = "IPCW", train.dat, 
	test.dat = NULL) { 
	
	quantile <- stats::quantile
	train.dat$eventsurv <- ifelse(train.dat$eventsurv.ind==TRUE, 1, 0)
	fmla <- as.formula(paste0(c(paste0("eventsurv~", paste0(varnms_fixed, collapse="+")),
		paste0(paste("s(", varnms_smooth[, 1],", by = ", varnms_smooth[, 2], ")", sep = ""), 
			collapse = "+")), collapse = "+"))
	
	if (cens_method == "IPCW") {
		GAM <- gam(fmla, data = train.dat, weights = wts, family = binomial)
	}	
	if (cens_method == "ZERO") {
		GAM <- gam(fmla, data = train.dat, family = binomial)
	}	
	if (cens_method == "DISCARD") {
		GAM <- gam(fmla, data = filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0)), 
			family = binomial)
	}
	
	 if(is.null(test.dat)) {
    probs <- predict(GAM, type = "response", na.action = na.omit)    
  } else {
    probs <- predict(GAM, newdata=test.dat, type = "response", na.action = na.omit)    
  }
  return(1-probs)
}





