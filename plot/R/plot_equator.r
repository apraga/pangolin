# Plot northern latitude closest to the equator to see how the error is
# decreasing 
library(ggplot2)

folder="/wkdir/pae2/praga/parallel/equator/"
vec <- c(160,320,640,800)
#vec <- c(20, 40, 80, 160,320,640)
q <- c()
lon <- c()
res <- c()
for (i in vec) {
  dat <- read.table(paste(folder, "extracted_", i, "_un_T.dat", sep=""))
  q <- c(q, dat$V1)
  lon <- c(lon, dat$V3)
  res <- c(res, rep(i, length(dat$V1)))
  i <- i*2
}

# Reference
dat <- read.table(paste(folder, "extracted_160_un_0.dat", sep=""))
q <- c(q, dat$V1)
lon <- c(lon, dat$V3)
res <- c(res, rep("ref", length(dat$V1)))

pdf("equator_cosine.pdf")
data <- data.frame(q, lon, res)
p <- ggplot(data, aes(x=lon, y=q, colour=factor(res))) + geom_line()
p <- p + scale_color_discrete(name="Nb lat/2")
p <- p + xlab(expression(lambda)) + ylab("ratio")
p <- p + ggtitle("Equator distribution for cosine bells")
print(p)
dev.off()
