criteria <- which(HP.test$SBP_Meds == 1)
criteria <- which(HP.test$SBP_Meds == 0)
criteria <- which(HP.test$SBP_Meds == 0 & HP.test$Has_Diab == 0)
criteria <- which(HP.test$SBP_Meds == 1 & HP.test$gender == 0 & HP.test$age <= 55)
criteria <- which(HP.test$age > 70)
criteria <- which(HP.test$SBP > 140)
criteria <- which(HP.test$SBP > 140 & HP.test$SBP_Meds == 1)
criteria <- which(HP.test$age > 68.5)
criteria <- which(HP.test$age >= 68.5 & HP.test$age< 79.5)
criteria <- which(HP.test$SBP_Meds == 1 & HP.test$gender == 0 & HP.test$age <= 55 & SBP)

criteria <- which(HP.test$SBP_Meds == 1 & HP.test$age <= 55)
criteria <- which(HP.test$age <= 55)
criteria <- which(HP.test$BMI > 25)
criteria <- which(HP.test$SBP_Meds == 1 & HP.test$age <= 55 & HP.test$Has_Diab == 0 )
criteria <- which(HP.test$SBP_Meds == 1 & HP.test$age <= 55 & HP.test$Has_Diab == 0 & HP.test$SBP < 140)
### pick this one
criteria <- which(HP.test$SBP_Meds == 1 & HP.test$age <= 55  & HP.test$SBP < 140)


ML.preds <- list(pRP[criteria], pLR[criteria], pGAM[criteria], pkNN[criteria], pBayes[criteria], pCPH[criteria])
risk.cutpts = c(Inf,1-c(0.025, 0.05, 0.075, 0.1),-Inf)


pred.rate <- 100*(1-unlist(lapply(ML.preds,mean)))
calib <- unlist(lapply(ML.preds,calib.stat,cutpts=risk.cutpts,test.dat=HP.test[criteria, ]))
c.index <- unlist(lapply(ML.preds,Cindex,test.dat=HP.test[criteria,]))

# Calibration and C-index Tables
rownms <- c("Tree", "Logistic", "GAM", "k-NN", "Bayes", "Cox")
colnms <- c("Predicted event rate (%)", "Calibration", "C-Index")
results.tab <- data.frame(pred.rate, calib, c.index)
rownames(results.tab) <- rownms
colnames(results.tab) <- colnms
results.tab

print( xtable(results.tab, align= "lccc",digits = c(0, 2, 2, 3)),
       hline.after=c(4, 8, 12, 16, 20) )
