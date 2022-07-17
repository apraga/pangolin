# Plot mixing diagnostics for 2 resolutions and with/without slope limiter

library(ggplot2)

folder <- "/wkdir/pae2/praga/parallel/mixing/"
n <- "40"
name <- paste(folder, "ratio_", n, "lat_sp_CFL0.7_0_cos.dat", sep="")
cos_init <- read.table(name)
name <- paste(folder, "ratio_", n, "lat_sp_CFL0.7_0_corr.dat", sep="")
corr_init <- read.table(name)

# Get limit
khi_min = min(cos_init$V1)
khi_max = max(cos_init$V1)
xi_min = min(corr_init$V1)
xi_max = max(corr_init$V1)


khi <- c()
xi <- c()
res <- c()
type <- c()
res <- c()
labels <- c("1.5°", "0.75°")
k <- 1

for (n in c("40", "80")) {
  for (case in c("sp", "un")) {
    name <- paste(folder, "ratio_", n, "lat_", case, "_CFL0.7_T2_cos.dat", sep="")
    cos <- read.table(name)
    name <- paste(folder, "ratio_", n, "lat_", case, "_CFL0.7_T2_corr.dat", sep="")
    corr <- read.table(name)
    khi <- c(khi, cos$V1)
    xi <- c(xi, corr$V1)
    res <- c(res, rep(labels[k], length(cos$V1)))
    type <- c(type, rep(case, length(cos$V1)))
  }
    k <- k +1
}

data <- data.frame(khi, xi, res, type)
# Change plot order
data$res <- factor(data$res, levels=c("1.5°", "0.75°"))

# Create correlation function
x <- seq(khi_min, khi_max, 0.001)
correl <- data.frame(x=x, y=-0.8*x**2+0.9)
# Create boundary (faster as a dataframe instead of manual geom_lines)
boundary <- data.frame(x=c(khi_min, khi_max, khi_max,khi_min), 
    y=c(xi_max, xi_max, xi_min, xi_max))

pdf("mixing.pdf", width=10,height=7)
# A bit of transparency helps
p <- ggplot(data, aes(x=khi, y=xi)) + geom_point(col="red", alpha=1/10)
p <- p + facet_wrap(type~res, ncol=2)  

p <- p + geom_line(data=correl, aes(x=x, y=y))
p <- p + geom_polygon(data=boundary, aes(x=x, y=y), col="black", alpha=0)
p <- p + coord_cartesian(xlim = c(0.05, 1.05), ylim = c(0,1)) 
p <- p + xlab(expression(chi)) + ylab(expression(xi))
p <- p + ggtitle("Mixing diagnostic")
print(p)
dev.off()
