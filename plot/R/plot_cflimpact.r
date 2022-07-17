# Plot cfl impact on accuracy
library(ggplot2)
library(reshape2)

dir <-"/wkdir/pae2/praga/parallel/cfl_impact/gaussian_160lat"
cur <- paste(dir, "/cfl-impact-un.dat", sep="")
err <- read.table(cur)
colnames(err) <- c("dt", "l2", "loo", "quad")

err2 <- melt(err, id='dt')
pdf("cfl_impact_gaussian_160lat.pdf")
p <- ggplot(err2, aes(x=dt, y=value, color=variable)) + geom_line() + geom_point()
p <- p + facet_grid(variable~.,scales="free")
#p <- p + scale_x_continuous(breaks=c(0.8,8,12,18,20,24))
p <- p + scale_x_continuous(breaks=c(0.6,3,5,7.5,9))
p <- p + xlab("dt (min)") + ylab("error")
p <- p + ggtitle("Timestep impact on accuracy (Gaussian, 2T, unlimited, 1.12x0.75°)")
print(p)
dev.off()
