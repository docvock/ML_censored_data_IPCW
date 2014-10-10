#  Bagged Recursive Partitioning Models
lossmat <- rbind(c(0,1),c(5,0))
set.seed(1101985)
pRP.bag <- rpart_bag_functions(varnms_a, 50, "IPCW", HP.train, HP.test, lossmat, "information",
                               cp = 0.00001, minbucket = 100)
pRP.bag.zero <- rpart_bag_functions(varnms_a, 50, "ZERO", HP.train, HP.test, lossmat, "information", 
                                    cp = 0.00001, minbucket = 100)
pRP.bag.disc <- rpart_functions(varnms_a, 50, "DISCARD", HP.train, HP.test, lossmat, "information", 
                                cp = 0.00001, minbucket = 100)

#  Bagged Logistic Regression Models
set.seed(1101985)
pLR.bag <- logistic_bag_functions(intx.varnms, 50, "IPCW", HP.train, HP.test)
pLR.bag.zero <- logistic_bag_functions(intx.varnms, 50, "ZERO", HP.train, HP.test)
pLR.bag.disc <- logistic_bag_functions(intx.varnms, 50, "DISCARD", HP.train, HP.test)

#  Bagged Generalized Additive Models
set.seed(1101985)
pGAM.bag <- GAM_bag_functions(varnms_fixed, varnms_smooth, 25, "IPCW", HP.train, 
                              HP.test)

#  Bagged MARS Models
set.seed(1101985)
pMARS <- MARS_bag_functions(varnms, 50, "IPCW", HP.train, HP.test)
pMARS.zero <- MARS_bag_functions(varnms, 50, "ZERO", HP.train, HP.test)
pMARS.disc <- MARS_bag_functions(varnms, 50, "DISC", HP.train, HP.test)

#  Bagged knn
set.seed(1101985)
pknn <- knn_bag_functions(varnms_stand, 50, "IPCW", HP.train, HP.test, k = 100)
pknn.zero <- knn_bag_functions(varnms_stand, 50, "ZERO", HP.train, HP.test, k = 100)
pknn.disc <- knn_bag_functions(varnms_stand, 50, "DISC", HP.train, HP.test, k = 100)