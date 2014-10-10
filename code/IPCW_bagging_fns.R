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
