KM.return <- function(t, surv.obj) {
    KM.return <- sapply(1:length(t), function(pos) {
    	surv.obj$surv[min(which(surv.obj$time > t[pos]))]
    })
  }

outer2 <- function(x) {
	outer2 <- outer(x,x)
	return(outer2)
}

## IPCW recursive partitioning

rpart_cvstats <- function(cp.seq, n.folds = 5, varnms, cens_method = "IPCW", losscost_seq, minbucket_seq, 
                          train.dat, split = "information", risk.cutpts = NULL, seed = 1,...) {  

  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
  
  set.seed(seed)
  allfolds <- createFolds(train.dat$eventsurv.ind, k = n.folds, list = FALSE)

  doFold <- function(j, lossmat, minbucket) {
    ind <- allfolds != j
    cvtrain.set <- subset(train.dat,ind)
    cvtest.set <- subset(train.dat,!ind)

    
    if(cens_method == "IPCW") {
      RP <- rpart(fmla, data = cvtrain.set, weights = wts, method="class", 
                  parms = list(loss = lossmat, split = split), 
      						control = rpart.control(cp=min(cp.seq), minbucket = minbucket, ...))
    }
    if(cens_method == "ZERO") {
      RP <- rpart(fmla, data = cvtrain.set, method = "class", 
                  parms = list(loss = lossmat, split = split), 
      						control = rpart.control(cp=min(cp.seq), minbucket = minbucket, ...))
    }
    if(cens_method == "DISCARD") {
      RP <- rpart(fmla, data = filter(cvtrain.set, !(T.use <= (SURVTIME*365) & C.use == 0)), 
                  method = "class", parms = list(loss = lossmat, split = split), 
      						control=rpart.control(cp=min(cp.seq), minbucket = minbucket, ...))
    } 
    
    if(cens_method == "SPLIT") { 
      sf <- survfit(Surv(T.use, C.use) ~ 1, data = cvtrain.set)
      sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
      cvtrain.set.nc <- filter(cvtrain.set,!(T.use <= (SURVTIME*365) & C.use == 0))
      cvtrain.set.nc <- mutate(cvtrain.set.nc, wts_split = 1) 	
      cvtrain.set.0 <- cvtrain.set.1 <- filter(cvtrain.set,(T.use <= (SURVTIME*365) & C.use == 0))
      cvtrain.set.0 <- mutate(cvtrain.set.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
      cvtrain.set.1 <- mutate(cvtrain.set.1, eventsurv.ind = TRUE, 
                          wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
      cvtrain.set.full <- rbind(cvtrain.set.nc, cvtrain.set.0, cvtrain.set.1)
      RP <- rpart(fmla, data = cvtrain.set.full, weights = wts_split, method = "class", 
                  parms = list(loss = lossmat, split = split), 
                  control = rpart.control(cp=min(cp.seq), minbucket = minbucket, ...))
    }
    
    cvstats <- t(sapply(cp.seq,function(mycp) {
      pruneRP <- prune(RP,mycp)
      preds <- predict(pruneRP, newdata=cvtest.set, type = "prob", na.action = na.omit)[,1]
      if(is.null(risk.cutpts)) cutpts <- quantile(preds) else { cutpts <- risk.cutpts }
      calib <- suppressWarnings(calib.stat(preds,cutpts=cutpts,test.dat=cvtest.set))
      cindex <- suppressWarnings(Cindex(preds,test.dat=cvtest.set))
      c(calib=calib,cindex=cindex,cp=mycp)
    }))
    
    return(data.frame(cvstats))
  }

   summ.table <- NULL
    
    for(i in 1:length(losscost_seq)) {
      for(k in 1:length(minbucket_seq)) {
        lossmat <- matrix(c(0, losscost_seq[i], 1, 0), nrow = 2, ncol = 2)
      	
  cvtable <- rbindlist(lapply(1:n.folds,doFold, lossmat = lossmat, minbucket = minbucket_seq[k]))
  tmp <- cvtable[,list(mean(calib),mean(cindex.concordant),sd(calib),sd(cindex.concordant)),by=cp]
  summ.table.ik <- data.frame(losscost = rep(losscost_seq[i], nrow(tmp)), 
                              minbucket = rep(minbucket_seq[k], nrow(tmp)), cp = tmp$cp,
                              calib = tmp$V1, cindex = tmp$V2, sd.calib = tmp$V3, sd.cindex = tmp$V4)
  summ.table <- rbind(summ.table, summ.table.ik)
      	print(c(i,k))
    } #close i loop
  } #close k loop

  return(summ.table)
  
}

rpart_functions <- function(varnms, cens_method = "IPCW", train.dat, test.dat = NULL, myseed = NULL, 
                            losscost_seq = c(2.5, 5, 10), minbucket_seq = c(50, 100, 200), split = "information", 
                            cp=0.00001, risk.cutpts = NULL, n.folds = 5, tuning.rule = "1SE",...) {  
  if(is.null(myseed) ==  TRUE) {
    myseed <- as.integer(Sys.time())
  }
  
  cp.seq <- 10^seq(-2,log10(cp),by=-0.05)
  
  cv.results <- rpart_cvstats(cp.seq, n.folds = n.folds, varnms = varnms, 
  	cens_method = cens_method, losscost_seq = losscost_seq, minbucket_seq = minbucket_seq, 
  	train.dat = train.dat, split = split, risk.cutpts = risk.cutpts, seed = myseed,...)
  
  ## Implement a modification of 1-SE rule pruning
  if(tuning.rule=="1SE") {
    #wmin.calib <- which.min(cv.results$calib)
    wmax.cindex <- which.max(cv.results$cindex)
    sdrule.cindex <- cv.results[which(cv.results$cindex >= (cv.results$cindex[wmax.cindex] - 
    		cv.results$sd.cindex[wmax.cindex])),]
    sdrule.cindex <- sdrule.cindex[which(sdrule.cindex$losscost == min(sdrule.cindex$losscost)), ]
  	sdrule.cindex <- sdrule.cindex[which(sdrule.cindex$minbucket == max(sdrule.cindex$minbucket)), ]


    final.cp <- sdrule.cindex$cp[which(sdrule.cindex$cp == max(sdrule.cindex$cp))]  ## Choose model
    final.minbucket <- sdrule.cindex$minbucket[which(sdrule.cindex$cp == max(sdrule.cindex$cp))]  ## Choose model
    final.losscost <- sdrule.cindex$losscost[which(sdrule.cindex$cp == max(sdrule.cindex$cp))]  ## Choose model
  }
  
  if(tuning.rule=="ranks") {
    rnks.calib <- rank(cv.results$calib)
    rnks.cindex <- rank(-cv.results$cindex)
    
    rnks <- (rnks.calib + rnks.cindex)/2
    min.rnk <- which.min(rnks)
    final.cp <- cv.results$cp[min.rnk]    
  }
  
  ##

  print("Cross-validation estimates of calibration and discrimination:")
  print(cv.results)
  print(paste0("Final CP: ",final.cp))
  print(paste0("Final Min Bucket: ",final.minbucket))
  print(paste0("Final Cost Ration: ",final.losscost))
  
  lossmat <- matrix(c(0, final.losscost, 1, 0), nrow = 2, ncol = 2)
  
  
	fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))      
  if(cens_method == "IPCW") {
    RP <- rpart(fmla, data = train.dat, weights = wts, method="class", 
                parms = list(loss = lossmat, split = split), 
    						control = rpart.control(cp=final.cp, minbucket = final.minbucket, ...))
  }
  if(cens_method == "ZERO") {
    RP <- rpart(fmla, data = train.dat, method = "class", 
                parms = list(loss = lossmat, split = split), 
    						control = rpart.control(cp=final.cp, minbucket = final.minbucket, ...))
  }
  if(cens_method == "DISCARD") {
    RP <- rpart(fmla, data = filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0)), 
                method = "class", parms = list(loss = lossmat, split = split), 
    						control=rpart.control(cp=final.cp, minbucket = final.minbucket, ...))
  }  
  if(cens_method == "SPLIT") {
    sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
    sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
    train.dat.nc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
    train.dat.nc <- mutate(train.dat.nc, wts_split = 1) 	
    train.dat.0 <- train.dat.1 <- filter(train.dat,(T.use <= (SURVTIME*365) & C.use == 0))
    train.dat.0 <- mutate(train.dat.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
    train.dat.1 <- mutate(train.dat.1, eventsurv.ind = TRUE, 
                          wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
    train.dat.full <- rbind(train.dat.nc, train.dat.0, train.dat.1)
    
    RP <- rpart(fmla, data = train.dat.full, weights = wts_split, method = "class", 
                parms = list(loss = lossmat, split = split), 
                control = rpart.control(cp=final.cp, minbucket = final.minbucket, ...))
  }

  if(is.null(test.dat)) {
    probs <- predict(RP, type = "prob", na.action = na.omit)[,1]    
  } else {
    probs <- predict(RP,newdata = test.dat, type = "prob", na.action = na.omit)[,1]
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
		paste0(paste("s(", varnms_smooth,")", sep = ""), 
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

## IPCW k-Nearest Neighbors

knn_cvstats <- function(neighbor.seq, n.folds = 5, varnms_stand, cens_method = "IPCW", dist_method = "mahalanobis",
                          train.dat, risk.cutpts = NULL, seed = 1) {  
  
  fmla <- as.formula(paste0("eventsurv.ind~",paste0(varnms,collapse="+")))
  
  set.seed(seed)
  allfolds <- createFolds(train.dat$eventsurv.ind, k = n.folds, list = FALSE)
  
  doFold <- function(j) {
    ind <- allfolds != j
    cvtrain.set <- subset(train.dat,ind)
    cvtest.set <- subset(train.dat,!ind)
    
    if (cens_method == "DISCARD") {
      cvtrain.set <- filter(cvtrain.set, !(T.use <= (SURVTIME*365) & C.use == 0))
    }
    if (cens_method == "SPLIT") {
      sf <- survfit(Surv(T.use, C.use) ~ 1, data = cvtrain.set)
      sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
      cvtrain.set <- mutate(cvtrain.set, wts_split = 1-sf.SURVTIME/KM.return(T.use, sf),
                          eventsurv = ifelse(T.use <= (SURVTIME*365) & C.use == 0, wts_split, eventsurv))
    }
    x <- rbind(cvtrain.set[, varnms_stand], cvtest.set[, varnms_stand])
    rownames(x) <- 1:nrow(x)
    y <- matrix(cvtrain.set$eventsurv, ncol = 1)
    y <- as.data.frame(y)
    colnames(y) <- c("Event")
    rownames(y) <- 1:nrow(cvtrain.set)
    
    msn <- yai(x = x,y = y, k = max(neighbor.seq), noRefs = TRUE, method = dist_method)
    if (cens_method == "IPCW") {
      neighbor.matrix <- matrix(as.numeric(msn$neiIdsTrgs),ncol = max(neighbor.seq), nrow = nrow(msn$neiIdsTrgs))
      knn.weight <- function(neighbor, max.neighbor, train.dat) {
        knn.weight <- weighted.mean(x = cvtrain.set$eventsurv[neighbor[1:max.neighbor]], 
                                    w = cvtrain.set$wts[neighbor[1:max.neighbor]])
        return(knn.weight)
      }
      probs <- NULL
      for(k in 1:length(neighbor.seq)) {
      probs <- cbind(probs, 1-apply(neighbor.matrix, 1, knn.weight, max.neighbor = neighbor.seq[k],
      	train.dat = cvtrain.set))
      }
    }
    if (cens_method == "ZERO" | cens_method == "DISCARD" | cens_method == "SPLIT" )  {
      neighbor.matrix <- matrix(as.numeric(msn$neiIdsTrgs),ncol = max(neighbor.seq), nrow = nrow(msn$neiIdsTrgs))
      knn.weight <- function(neighbor, max.neighbor, train.dat) {
        knn.weight <- mean(x = cvtrain.set$eventsurv[neighbor[1:max.neighbor]])
        return(knn.weight)
      }
      probs <- NULL
      for(k in 1:length(neighbor.seq)) {
        probs <- cbind(probs, 1-apply(neighbor.matrix, 1, knn.weight, max.neighbor = neighbor.seq[k],
      		train.dat = cvtrain.set))
      }
    }  	
      
        cvstats <- apply(probs, 2, function(preds) {
          if(is.null(risk.cutpts)) cutpts <- quantile(preds) else { cutpts <- risk.cutpts }
          calib <- suppressWarnings(calib.stat(preds, cutpts = cutpts, test.dat = cvtest.set))
          cindex <- suppressWarnings(Cindex(preds, test.dat = cvtest.set))
          c(calib = calib, cindex = cindex)
        })
        
    cvstats <- as.data.frame(t(cvstats))
    cvstats$neighbor <- neighbor.seq
        return(cvstats)
      }
      
      cvtable <- rbindlist(lapply(1:n.folds,doFold))
      tmp <- cvtable[,list(mean(calib),mean(cindex.concordant),sd(calib),sd(cindex.concordant)),by=neighbor]
      summ.table <- data.frame(neighbor = tmp$neighbor, calib = tmp$V1, cindex = tmp$V2, 
      	sd.calib = tmp$V3, sd.cindex = tmp$V4)
  
  return(summ.table)
  
}

knn_functions <- function(varnms_stand, kmax = 500, cens_method = "IPCW", dist_method = "mahalanobis", 
                          train.dat, test.dat, risk.cutpts = NULL, myseed = NULL, n.folds = 5, 
													tuning.rule = "1SE") { 
  train.dat$eventsurv <- ifelse(train.dat$eventsurv.ind==TRUE, 1, 0)
  neighbor.seq <- seq(from = 50, to = kmax, by = 50)
  
  cv.results <- knn_cvstats(neighbor.seq = neighbor.seq, n.folds = 5, 
  	varnms_stand = varnms_stand, cens_method = cens_method, dist_method = dist_method,
  	train.dat = train.dat, risk.cutpts = risk.cutpts, seed = myseed)
  if(tuning.rule=="1SE") {
    wmax.cindex <- which.max(cv.results$cindex)
    sdrule.cindex <- cv.results[which(cv.results$cindex >= (cv.results$cindex[wmax.cindex] - 
                                                              cv.results$sd.cindex[wmax.cindex])),]
    
    final.neighbor <- sdrule.cindex$neighbor[which(sdrule.cindex$neighbor == 
    		max(sdrule.cindex$neighbor))]  ## Choose model
  }
  
  print("Cross-validation estimates of calibration and discrimination:")
  print(cv.results)
  print(paste0("Final Neighbor: ",final.neighbor))
  
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
  
  msn <- yai(x = x,y = y, k = final.neighbor, noRefs = TRUE, method = dist_method)
  if (cens_method == "IPCW") {
    neighbor.matrix <- matrix(as.numeric(msn$neiIdsTrgs),ncol = final.neighbor, nrow = nrow(msn$neiIdsTrgs))
    knn.weight <- function(neighbor, train.dat) {
      knn.weight <- weighted.mean(x = train.dat$eventsurv[neighbor], w = train.dat$wts[neighbor])
      return(knn.weight)
    }
    probs <- apply(neighbor.matrix, 1, knn.weight, train.dat)
    return(1-probs)
  }
  if (cens_method == "ZERO" | cens_method == "DISCARD" | cens_method == "SPLIT" )  {
  	  neighbor.matrix <- matrix(as.numeric(msn$neiIdsTrgs),ncol = final.neighbor, nrow = nrow(msn$neiIdsTrgs))
      knn.weight <- function(neighbor, train.dat) {
        knn.weight <- mean(x = train.dat$eventsurv[neighbor])
        return(knn.weight)
      }
    probs <- apply(neighbor.matrix, 1, knn.weight, train.dat)
    return(1-probs)
  }		
}

## IPCW Bayesian Network Functions

Bayes_functions <- function(cens_method = "IPCW", train.dat, test.dat = NULL, 
	ID = 1:nrow(train.dat),...) {  
	
	if(cens_method == "IPCW") {
		IPCW_design <- svydesign(id = ~ (ID), weights = ~wts, data = train.dat)
	}
	if(cens_method == "ZERO") {
		IPCW_design <- svydesign(id = ~ (ID), data = train.dat)
	}
	if(cens_method == "DISCARD") {
		train.dat.disc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		ID <- 1:nrow(train.dat.disc)
		IPCW_design <- svydesign(id = ~ (ID), data = train.dat.disc)
	}	
	if(cens_method == "SPLIT") {
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		sf.SURVTIME <- sf$surv[min(which(sf$time > SURVTIME*365))]
		train.dat.nc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.nc <- mutate(train.dat.nc, wts_split = 1) 	
		train.dat.0 <- train.dat.1 <- filter(train.dat,(T.use <= (SURVTIME*365) & C.use == 0))
		train.dat.0 <- mutate(train.dat.0, wts_split = sf.SURVTIME/KM.return(T.use, sf))
		train.dat.1 <- mutate(train.dat.1, eventsurv.ind = TRUE, 
		wts_split = 1-sf.SURVTIME/KM.return(T.use,sf))
		train.dat.full <- rbind(train.dat.nc, train.dat.0, train.dat.1)
		
		IPCW_design <- svydesign(id = ~ (ID), weights = ~wts_split, data = train.dat.full)
	}
	
	p.age_event <- prop.table(svytable(~ eventsurv.ind + age_f, IPCW_design), 1)
	p.BMI_age.event <- prop.table(svytable(~ eventsurv.ind + age_f + BMI_f, IPCW_design), c(1, 2))
	p.diab_BMI.event <- prop.table(svytable(~ eventsurv.ind + BMI_f + Has_Diab_f, IPCW_design), c(1, 2))
	p.SBPMeds_age.BMI.diab.event <- prop.table(svytable(~ eventsurv.ind + age_f + BMI_f + Has_Diab_f + 
		SBP_Meds_f, IPCW_design), c(1, 2, 3, 4))
	p.smoke_age.event <- prop.table(svytable(~ eventsurv.ind + age_f + Smoking_f, IPCW_design), c(1, 2))

	if(cens_method == "IPCW") {
		p.SBP_age.BMI.smoke.SBPMeds.event <- lm(log_SBP ~ eventsurv.ind + SBP_Meds_f + 
				SBP_Meds_f*as.factor(BMI_f) + as.factor(age_f) + as.factor(Smoking), weights = wts, 
				data = train.dat)
		sigma_SBP <- sqrt(weighted.mean(resid(p.SBP_age.BMI.smoke.SBPMeds.event)^2, train.dat$wts))
		p.HDL_age.BMI.diab.event <- lm(log_HDL ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab  , weights = wts, data = train.dat)
		p.TC_age.BMI.diab.event <- lm(log_TC ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab  , weights = wts, data = train.dat)
		resid.HDL.TC <- cbind(resid(p.HDL_age.BMI.diab.event),resid(p.TC_age.BMI.diab.event))
		vcov.resid.HDL.TC <- t(apply(resid.HDL.TC, 1, outer2))
		vcov.resid.HDL.TC <- apply(vcov.resid.HDL.TC, 2, weighted.mean, train.dat$wts)
		vcov.resid.HDL.TC <- matrix(vcov.resid.HDL.TC, nrow = 2, ncol = 2) 
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		p.noevent <- sf$surv[min(which(sf$time > SURVTIME*365))]
		p.event <- 1 - p.noevent
	}
	if(cens_method == "ZERO") {
		p.SBP_age.BMI.smoke.SBPMeds.event <- lm(log_SBP ~ eventsurv.ind + SBP_Meds_f + 
				SBP_Meds_f*as.factor(BMI_f) + as.factor(age_f) + as.factor(Smoking),  
				data = train.dat)
		sigma_SBP <- sqrt(mean(resid(p.SBP_age.BMI.smoke.SBPMeds.event)^2))
		p.HDL_age.BMI.diab.event <- lm(log_HDL ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab,  data = train.dat)
		p.TC_age.BMI.diab.event <- lm(log_TC ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab,  data = train.dat)
		resid.HDL.TC <- cbind(resid(p.HDL_age.BMI.diab.event), resid(p.TC_age.BMI.diab.event))
		vcov.resid.HDL.TC <- t(apply(resid.HDL.TC, 1, outer2))
		vcov.resid.HDL.TC <- apply(vcov.resid.HDL.TC, 2, mean)
		vcov.resid.HDL.TC <- matrix(vcov.resid.HDL.TC, nrow = 2, ncol = 2) 
		p.noevent <- mean(train.dat$eventsurv.ind == FALSE)
		p.event <- 1 - p.noevent
	}
	if(cens_method == "DISCARD") {
		train.dat.disc <- filter(train.dat,!(T.use <= (SURVTIME*365) & C.use == 0))
		p.SBP_age.BMI.smoke.SBPMeds.event <- lm(log_SBP ~ eventsurv.ind + SBP_Meds_f + 
				SBP_Meds_f*as.factor(BMI_f) + as.factor(age_f) + as.factor(Smoking),  
				data = train.dat.disc )
		sigma_SBP <- sqrt(mean(resid(p.SBP_age.BMI.smoke.SBPMeds.event)^2))
		p.HDL_age.BMI.diab.event <- lm(log_HDL ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab,  data = train.dat.disc )
		p.TC_age.BMI.diab.event <- lm(log_TC ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab,  data = train.dat.disc )
		resid.HDL.TC <- cbind(resid(p.HDL_age.BMI.diab.event), resid(p.TC_age.BMI.diab.event))
		vcov.resid.HDL.TC <- t(apply(resid.HDL.TC, 1, outer2))
		vcov.resid.HDL.TC <- apply(vcov.resid.HDL.TC, 2, mean)
		vcov.resid.HDL.TC <- matrix(vcov.resid.HDL.TC, nrow = 2, ncol = 2) 
		p.noevent <- mean(train.dat.disc$eventsurv.ind == FALSE)
		p.event <- 1 - p.noevent
	}
	if(cens_method == "SPLIT") {
		p.SBP_age.BMI.smoke.SBPMeds.event <- lm(log_SBP ~ eventsurv.ind + SBP_Meds_f + 
				SBP_Meds_f*as.factor(BMI_f) + as.factor(age_f) + as.factor(Smoking), weights = wts_split, 
				data = train.dat.full)
		sigma_SBP <- sqrt(weighted.mean(resid(p.SBP_age.BMI.smoke.SBPMeds.event)^2, train.dat.full$wts_split))
		p.HDL_age.BMI.diab.event <- lm(log_HDL ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab  , weights = wts_split, data = train.dat.full)
		p.TC_age.BMI.diab.event <- lm(log_TC ~ eventsurv.ind + as.factor(BMI_f) + as.factor(age_f) + 
				Has_Diab  , weights = wts_split, data = train.dat.full)
		resid.HDL.TC <- cbind(resid(p.HDL_age.BMI.diab.event),resid(p.TC_age.BMI.diab.event))
		vcov.resid.HDL.TC <- t(apply(resid.HDL.TC, 1, outer2))
		vcov.resid.HDL.TC <- apply(vcov.resid.HDL.TC, 2, weighted.mean, train.dat.full$wts_split)
		vcov.resid.HDL.TC <- matrix(vcov.resid.HDL.TC, nrow = 2, ncol = 2) 
		sf <- survfit(Surv(T.use, C.use) ~ 1, data = train.dat)
		p.noevent <- sf$surv[min(which(sf$time > SURVTIME*365))]
		p.event <- 1 - p.noevent
	}
	
	if(is.null(test.dat)) {
		pred.dat <- train.dat 
	} else {
		pred.dat <- test.dat
	}
	
	pred.dat.1 <- pred.dat.0 <- pred.dat
	pred.dat.1$eventsurv.ind <- as.factor(rep(TRUE, nrow(pred.dat))) 
	pred.dat.0$eventsurv.ind <- as.factor(rep(FALSE, nrow(pred.dat)))
	pred.HDL.1 <- predict(p.HDL_age.BMI.diab.event, newdata=pred.dat.1)
	pred.HDL.0 <- predict(p.HDL_age.BMI.diab.event, newdata=pred.dat.0)

	pred.TC.1 <- predict(p.TC_age.BMI.diab.event, newdata=pred.dat.1)
	pred.TC.0 <- predict(p.TC_age.BMI.diab.event, newdata=pred.dat.0)

p.X.event <- 	p.age_event[2, pred.dat$age_f] *
	p.BMI_age.event[cbind(2, pred.dat$age_f, pred.dat$BMI_f)] *	
	p.diab_BMI.event[cbind(2, pred.dat$BMI_f, pred.dat$Has_Diab_f)] *
	p.SBPMeds_age.BMI.diab.event[cbind(2, pred.dat$age_f, pred.dat$BMI_f, pred.dat$Has_Diab_f, 
		pred.dat$SBP_Meds_f)] *
	p.smoke_age.event[cbind(2, pred.dat$age_f, pred.dat$Smoking_f)] *  
	dnorm(pred.dat$log_SBP, mean = predict(p.SBP_age.BMI.smoke.SBPMeds.event, newdata=pred.dat.1),
		sd = sigma_SBP) * 
	sapply(1:nrow(pred.dat), function(i) {
		dmvnorm(c(pred.dat$log_HDL[i], pred.dat$log_TC[i]), 
		c(pred.HDL.1[i], pred.TC.1[i]), sigma = vcov.resid.HDL.TC) })

p.X.noevent <- 	p.age_event[1, pred.dat$age_f] *
	p.BMI_age.event[cbind(1, pred.dat$age_f, pred.dat$BMI_f)] *	
	p.diab_BMI.event[cbind(1, pred.dat$BMI_f, pred.dat$Has_Diab_f)] *
	p.SBPMeds_age.BMI.diab.event[cbind(1, pred.dat$age_f, pred.dat$BMI_f, pred.dat$Has_Diab_f, 
		pred.dat$SBP_Meds_f)] *
	p.smoke_age.event[cbind(1, pred.dat$age_f, pred.dat$Smoking_f)] *  
	dnorm(pred.dat$log_SBP, mean = predict(p.SBP_age.BMI.smoke.SBPMeds.event, newdata=pred.dat.0),
		sd = sigma_SBP) * 
	sapply(1:nrow(pred.dat), function(i) {
		dmvnorm(c(pred.dat$log_HDL[i], pred.dat$log_TC[i]), 
		c(pred.HDL.0[i], pred.TC.0[i]), sigma = vcov.resid.HDL.TC) })

pBayes <- p.X.noevent*p.noevent/(p.X.noevent*p.noevent+p.X.event*p.event)
return(pBayes)
}

