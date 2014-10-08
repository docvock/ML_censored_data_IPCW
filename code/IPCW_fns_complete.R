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

# Recursive Partitioning Bagging Functions
rpart_bag_functions <- function(varnms, B, cens_method = "IPCW", train.dat.all, test.dat = NULL,
	lossmat = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), split = "information",...) { 

	if(is.null(test.dat)) {
		pred.prob <- matrix(0, nrow = nrow(train.dat.all), ncol = B)
	} else {
		pred.prob <- matrix(0, nrow = nrow(test.dat), ncol = B)	
	}
	fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))

	for(q in 1:B) {
		if(cens_method == "IPCW") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE, weight = wts)
		}
		if(cens_method == "ZERO") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE)
		}
		if(cens_method == "DISCARD") {
			train.dat.disc <- filter(train.dat.all, !(T.use <= (SURVTIME*365) & C.use == 0))
			train.dat <- sample_n(train.dat.disc, nrow(train.dat.all), replace=TRUE)
		}
		RP_B <- rpart(fmla, data = train.dat, method = "class", 
			parms = list(loss = lossmat, split = split), control = rpart.control(...))
		
		if(is.null(test.dat)) {
			pred.prob.B <- predict(RP_B, type = "prob", na.action=na.omit)[,1]    
		} else {
			pred.prob.B <- predict(RP_B, newdata=test.dat, type = "prob", na.action=na.omit)[,1]    
  	}
		pred.prob[, q] <- pred.prob.B
		print(q)
	}
	probs <- apply(pred.prob, 1, mean)
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
    probs <- as.vector(predict(LR, type = "response", na.action=na.omit)    
  } else {
    probs <- predict(LR, newdata=test.dat, type = "response", na.action = na.omit)
  }
  return(1-probs)
}

# Logistic Bagging Functions
logistic_bag_functions <- function(varnms, B, cens_method = "IPCW", train.dat.all, test.dat = NULL) { 

	if(is.null(test.dat)) {
		pred.prob <- matrix(0, nrow = nrow(train.dat.all), ncol = B)
	} else {
		pred.prob <- matrix(0, nrow = nrow(test.dat), ncol = B)	
	}

  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))

	for(q in 1:B) {
		if(cens_method == "IPCW") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE, weight = wts)
		}
		if(cens_method == "ZERO") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE)
		}
		if(cens_method == "DISCARD") {
			train.dat.disc <- filter(train.dat.all, !(T.use <= (SURVTIME*365) & C.use == 0))
			train.dat <- sample_n(train.dat.disc, nrow(train.dat.all), replace=TRUE)
		}
		LR_B <- glm(fmla, data = train.dat,family = "binomial")
		
		if(is.null(test.dat)) {
			pred.prob.B <- predict(LR_B, type = "response", na.action=na.omit)    
		} else {
			pred.prob.B <- predict(LR_B, newdata=test.dat, type = "response", na.action=na.omit)    
  	}
		pred.prob[, q] <- pred.prob.B
		print(q)
	}
	probs <- apply(pred.prob, 1, mean)
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

# Generalized Additive Model Bagging Functions
GAM_bag_functions <- function(varnms_fixed, varnms_smooth, B, cens_method = "IPCW", train.dat.all,
	test.dat = NULL) { 

	if(is.null(test.dat)) {
		pred.prob <- matrix(0, nrow = nrow(train.dat.all), ncol = B)
	} else {
		pred.prob <- matrix(0, nrow = nrow(test.dat), ncol = B)	
	}
	quantile <- stats::quantile
	train.dat.all$eventsurv <- ifelse(train.dat.all$eventsurv.ind==TRUE, 1, 0)
	fmla <- as.formula(paste0(c(paste0("eventsurv~", paste0(varnms_fixed, collapse="+")),
		paste0(paste("s(", varnms_smooth[, 1],", by = ", varnms_smooth[, 2], ")", sep = ""), 
			collapse = "+")), collapse = "+"))

	for(q in 1:B) {
		if(cens_method == "IPCW") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE, weight = wts)
		}
		if(cens_method == "ZERO") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE)
		}
		if(cens_method == "DISCARD") {
			train.dat.disc <- filter(train.dat.all, !(T.use <= (SURVTIME*365) & C.use == 0))
			train.dat <- sample_n(train.dat.disc, nrow(train.dat.all), replace=TRUE)
		}
		GAM_B <- gam(fmla, data = train.dat, family = binomial)
		
		if(is.null(test.dat)) {
			pred.prob.B <- predict(GAM_B, type = "response", na.action=na.omit)    
		} else {
			pred.prob.B <- predict(GAM_B, newdata=test.dat, type = "response", na.action=na.omit)    
  	}
		pred.prob[, q] <- pred.prob.B
		print(q)
	}
	probs <- apply(pred.prob, 1, mean)
	return(1-probs)
	}

# MARS Bagging Functions
MARS_bag_functions <- function(varnms, B, cens_method = "IPCW", train.dat.all, test.dat = NULL) { 

	if(is.null(test.dat)) {
		pred.prob <- matrix(0, nrow = nrow(train.dat.all), ncol = B)
	} else {
		pred.prob <- matrix(0, nrow = nrow(test.dat), ncol = B)	
	}
	train.dat.all$eventsurv <- ifelse(train.dat.all$eventsurv.ind==TRUE, 1, 0)
  fmla <- as.formula(paste0("eventsurv~",paste0(varnms,collapse="+")))

	for(q in 1:B) {
		if(cens_method == "IPCW") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE, weight = wts)
		}
		if(cens_method == "ZERO") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE)
		}
		if(cens_method == "DISCARD") {
			train.dat.disc <- filter(train.dat.all, !(T.use <= (SURVTIME*365) & C.use == 0))
			train.dat <- sample_n(train.dat.disc, nrow(train.dat.all), replace=TRUE)
		}
		MARS_B <- earth(fmla, glm=list(family=binomial), nk = 50, degree = 3, thresh = 0.0005, 
				data=train.dat)
		if(is.null(test.dat)) {
			pred.prob.B <- predict(MARS_B, type = "response")    
		} else {
			pred.prob.B <- predict(MARS_B, newdata=test.dat, type = "response")    
  	}
		pred.prob[, q] <- pred.prob.B
		print(q)
	}
	probs <- apply(pred.prob, 1, mean)
	return(1-probs)
	}

# kNN Bagging Functions
knn_bag_functions <- function(varnms, B, cens_method = "IPCW", train.dat.all, test.dat = NULL,...) { 

	if(is.null(test.dat)) {
		pred.prob <- matrix(0, nrow = nrow(train.dat.all), ncol = B)
	} else {
		pred.prob <- matrix(0, nrow = nrow(test.dat), ncol = B)	
	}
  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))

	for(q in 1:B) {
		if(cens_method == "IPCW") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE, weight = wts)
		}
		if(cens_method == "ZERO") {
			train.dat <- sample_n(train.dat.all, nrow(train.dat.all), replace=TRUE)
		}
		if(cens_method == "DISCARD") {
			train.dat.disc <- filter(train.dat.all, !(T.use <= (SURVTIME*365) & C.use == 0))
			train.dat <- sample_n(train.dat.disc, nrow(train.dat.all), replace=TRUE)
		}
		knn_B <- knn3(fmla,  data=train.dat,...)
		if(is.null(test.dat)) {
			pred.prob.B <- predict(knn_B, type = "prob")[, 2]    
		} else {
			pred.prob.B <- predict(knn_B, newdata=test.dat, type = "prob")[, 2]    
  	}
		pred.prob[, q] <- pred.prob.B
		print(q)
	}
	probs <- apply(pred.prob, 1, mean)
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
	test.dat$p <- 1-p
	Cindex_all <- survConcordance(Surv(T.use, C.use) ~ p, data = test.dat)
	Cindex <- Cindex_all$concordance
	return(Cindex)
}
	

