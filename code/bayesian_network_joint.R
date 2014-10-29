library(igraph)
library(shape)

Edges <- data.frame(
    from = c(rep("Event", 8), rep("Gender", 3), rep("Age", 5), rep("BMI", 4), rep("Diab", 2),
    	rep("SBP Med", 1), rep("Smoke", 1)),
    to = c("Gender", "Age", "BMI", "Diab", "SBP Med", "Smoke", "SBP", "HDL, TC", 
    				"Age", "SBP", "HDL, TC", 
    				"BMI", "SBP Med", "Smoke", "SBP", "HDL, TC",
    				"Diab", "SBP Med", "SBP", "HDL, TC",
    				"SBP Med", "HDL, TC",
    				"SBP",
    				"SBP"
    				))

g <- graph.edgelist(as.matrix(Edges))
coord.use <- as.matrix(read.table("coord_use.txt", sep = "\t", header = F))
coord.use[1, 1] <- (min(coord.use[, 1]) + max(coord.use[, 2]))/2  - 20 
coord.use[7,1] <- coord.use[7,1] + 5
coord.use <- coord.use[-9, ]


curv.use <- rep(0, 24)
curv.use[10] <- -1.8
curv.use[16] <- 2.1
vertex.size.use <- rep(20, 9)+5
vertex.size.use[6] <- 40
vertex.size.use[2] <- 30+5
vertex.size.use[7] <- 25+7
vertex.size.use[1] <- 25+5
vertex.size.use[9] <- 40

par(mar=c(2,0.1,0,0.1))
plot(g, vertex.shape="rectangle", vertex.color="white", vertex.size = vertex.size.use, 
	vertex.label.color = "black", edge.color = c(rep("white",8), rep("black", 16)), layout = coord.use, 
	margin = c(0, 0, 0, 0), 
	edge.curved = curv.use, 
	vertex.label.cex = 0.9, margin = c(0, 0, 0, 0))
rect(-1.25, -1.45, 1.2, 0.6, lwd = 2)
Arrows(-0.03, 0.9, -0.03, 0.67, lwd = 1, arr.length = 0.6)


id <- tkplot(g, canvas.width=450, canvas.height=500, vertex.shape="rectangle", vertex.color="white")
t4 <- tkplot.getcoords(id)
write.table(t4, "coord_use.txt", sep="\t", row.names = F, col.names = F)

