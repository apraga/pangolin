# Plot true and approximate nb cells on a latitude
library(ggplot2)
library(reshape2) # for melt
library(grid) # for unit
lats <- seq(89.5,0.5)
y_approx <- c()
y_true <- c()
y_trunc <- c()
coef <- pi/180.
dlat <- 1*coef

n <- length(lats)
for (i in seq( n)) {
  y_approx[i] <- 3*(2*i-1)
  y_trunc[i] <- min(3*(2*i-1), 344)
  y_true[i] <- 3*(2*sin(0.5*coef)*sin((i-0.5)*coef))/(1-cos(coef))
}

#postscript("nbcells.ps",width=10,height=10,pointsize=15)
#type <- c(rep(1,n), rep(2,n))
data <- data.frame(lats, y_approx, y_trunc, y_true)
data2 <- melt(data, id="lats")
levels(data2$variable) <- c("Linear", "Truncated", "Exact")
p <- ggplot(data2, aes(x=lats,y=value)) + geom_line(aes(linetype=variable),
                                                        size=1.5)
p <- p + xlab("Latitude") + ylab("Nb cells")
# Larger text, no legend title
p <- p + theme(text = element_text(size=28)) 
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1),
               legend.title=element_blank(),
               legend.key.width=unit(1,"cm"),
               legend.key.height=unit(1,"cm"))
               
# Reverse axis
p <- p + scale_x_reverse()
p <- p + ggtitle("Nb cells approximation (1x0.67°)")
ggsave("nb_cells.pdf", width=8, height=8)
#print(p)
