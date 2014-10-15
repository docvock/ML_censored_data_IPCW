stand <- function(x){
	stand <- (x-mean(x))/sd(x)
}

KM.return <- function(t, surv.obj) {
    KM.return <- sapply(1:length(t), function(pos) {
    	surv.obj$surv[min(which(surv.obj$time > t[pos]))]
    })
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
	if (cens_method == "SPLIT") {
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
		train.dat.nc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.nc <- mutate(train.dat.nc, wts_split = 1) 	
		train.dat.0 <- train.dat.1 <- filter(train.dat,(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.0 <- mutate(train.dat.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
		train.dat.1 <- mutate(train.dat.1, eventsurv.ind = TRUE, 
			wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
		train.dat.full <- rbind(train.dat.nc, train.dat.0, train.dat.1)

		RP <- rpart(fmla, data = train.dat.full, weights = wts_split, method="class", 
			parms = list(loss = lossmat, split = split), control=rpart.control(...))
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
	if (cens_method == "SPLIT") {
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
		train.dat.nc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.nc <- mutate(train.dat.nc, wts_split = 1) 	
		train.dat.0 <- train.dat.1 <- filter(train.dat,(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.0 <- mutate(train.dat.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
		train.dat.1 <- mutate(train.dat.1, eventsurv.ind = TRUE, 
			wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
		train.dat.full <- rbind(train.dat.nc, train.dat.0, train.dat.1)
		
		LR <- glm(fmla, data = train.dat.full, weight = wts_split, family = "binomial")
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
	if (cens_method == "SPLIT") {
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
		train.dat.nc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.nc <- mutate(train.dat.nc, wts_split = 1) 	
		train.dat.0 <- train.dat.1 <- filter(train.dat,(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.0 <- mutate(train.dat.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
		train.dat.1 <- mutate(train.dat.1, eventsurv = 1, 
			wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
		train.dat.full <- rbind(train.dat.nc, train.dat.0, train.dat.1)
		
		GAM <- gam(fmla, data = train.dat.full, weights = wts_split, family = binomial)
	}	
	
	 if(is.null(test.dat)) {
    probs <- predict(GAM, type = "response", na.action = na.omit)    
  } else {
    probs <- predict(GAM, newdata=test.dat, type = "response", na.action = na.omit)    
  }
  return(1-probs)
}



# IPCW knn
knn_functions <- function(varnms_stand, k = 100, cens_method = "IPCW", dist_method = "mahalanobis", 
	train.dat, test.dat) { 
	train.dat$eventsurv <- ifelse(train.dat$eventsurv.ind==TRUE, 1, 0)
	if (cens_method == "DISCARD") {
		train.dat <- filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0))
	}
	if (cens_method == "SPLIT") {
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
		train.dat <- mutate(train.dat, wts_split = 1-sf.SURVTIME/KM.return(T.use, sf),
			eventsurv = ifelse(T.use <= (SURVTIME*365) & C.use == 0, wts_split, eventsurv))
	}
	x <- rbind(train.dat[, varnms_stand], test.dat[, varnms_stand])
	rownames(x) <- 1:nrow(x)
	y <- matrix(train.dat$eventsurv, ncol = 1)
	y <- as.data.frame(y)
	colnames(y) <- c("Event")
	rownames(y) <- 1:nrow(train.dat)
	
	msn <- yai(x = x,y = y, k = k, noRefs = TRUE, method = dist_method)
	if (cens_method == "IPCW") {
		neighbor.matrix <- matrix(as.numeric(msn$neiIdsTrgs),ncol=100, nrow = nrow(msn$neiIdsTrgs))
		knn.weight <- function(neighbor, train.dat) {
			knn.weight <- weighted.mean(x = train.dat$eventsurv[neighbor], w = train.dat$wts[neighbor])
			return(knn.weight)
		}
		probs <- apply(neighbor.matrix, 1, knn.weight, train.dat)
		return(1-probs)
	}
	if (cens_method == "ZERO" | cens_method == "DISCARD" | cens_method == "SPLIT" )  {
		probs <- impute(msn, method = "mean")
		return(1-probs[, 1])
	}		
	
}

