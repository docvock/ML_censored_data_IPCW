library(igraph)
library(shape)

Edges <- data.frame(
    from = c(rep("Event", 9), rep("Gender", 4), rep("Age", 6), rep("BMI", 5), rep("Diab", 3),
    	rep("SBP Med", 1), rep("Smoke", 1)),
    to = c("Gender", "Age", "BMI", "Diab", "SBP Med", "Smoke", "SBP", "HDL", "TC",
    				"Age", "SBP", "HDL", "TC",
    				"BMI", "SBP Med", "Smoke", "SBP", "HDL", "TC",
    				"Diab", "SBP Med", "SBP", "HDL", "TC",
    				"SBP Med", "HDL", "TC",
    				"SBP",
    				"SBP"
    				))

g <- graph.edgelist(as.matrix(Edges))
coord.use <- as.matrix(read.table("coord_use.txt", sep = "\t", header = F))
coord.use[7,1] <- coord.use[7,1] + 5
curv.use <- rep(0, 29)
curv.use[12] <- -2.5
curv.use[11] <- -2
curv.use[19] <- 2.1
curv.use[18] <- -0.5
vertex.size.use <- rep(20, 10)+5
vertex.size.use[6] <- 40+10
vertex.size.use[2] <- 30+10
vertex.size.use[7] <- 25+7
vertex.size.use[1] <- 25+5

#vertex.size.use <- rep(20, 10)
#vertex.size.use[6] <- 40
#vertex.size.use[2] <- 30
#vertex.size.use[7] <- 25
#vertex.size.use[1] <- 25



plot(g, vertex.shape="rectangle", vertex.color="white", vertex.size = vertex.size.use, 
	vertex.label.color = "black", edge.color = c(rep("white",9), rep("black", 20)), layout = coord.use, 
	margin = c(0, 0, 0, 0), 
	edge.curved = curv.use, 
	vertex.label.cex = 0.9)
rect(-1.6, -1.6, 1.2, 0.6, lwd = 2)
Arrows(-0.24, 0.9, -0.24, 0.67, lwd = 1, arr.length = 0.6)


id <- tkplot(g, canvas.width=450, canvas.height=500, vertex.shape="rectangle", vertex.color="white")
t4 <- tkplot.getcoords(id)
write.table(t4, "coord_use.txt", sep="\t", row.names = F, col.names = F)

