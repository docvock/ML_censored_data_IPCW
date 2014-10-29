
hist(HP$DaysToEvent_Fram/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects in Cohort", main = "", col = "grey")
hist(HP$DaysToEvent_Fram[HP$CVDEvent_Fram == 1]/365, breaks = seq(from = 0, to =10, by = 0.5), 
	xlab = "Years from Index Date", ylab = "Number of Subjects in Cohort", main = "", col = "white", add = TRUE)

t <- seq(from = 0, to = 8, by = 0.5)
unknown.prop <- NULL
for (i in 1:length(t)) {
unknown <- 1-prop.table(table(ifelse(HP$T.use < (t[i]*365) & HP$C.use==0, 1, 0)))[1]
unknown.prop <- c(unknown.prop, unknown)
	print(i)
}

plot(t, unknown.prop, type = "b", pch = 19, ylab = "Proportion with Event Status Unknown", 
	xlab = "Years from Index Date")



dim(HP)
table(HP$gender)
