# Plot mixing diagnostics for 2 resolutions and with/without slope limiter

library(ggplot2)

source("mixing.r")

folder_pango <- "/wkdir/pae2/praga/parallel/mixing/"
folder <- "/wkdir/pae2/praga/lauritzen/scatter/data/"

names <- c( "mix-SLFV-ML-dx0.75-sp-CFL0.8.dat",
           "mix-UCISOM-dx0.75-sp-CFL1.0.dat",
           "mix-CAM-FV-dx0.75-sp-CFL1.2.dat")

models <- c("SLFV-ML CFL0.8", "UCISOM CFL1.0", "CAM-FV CFL1.2")

khi <- c()
xi <- c()
model <- c()
type <- c()

create_data()

data <- data.frame(khi, xi, model, type)
data$type <- factor(data$type, levels=c("unlimited", "shape preserving"))
# Create correlation function
khi_min <- 0.1
khi_max <- 1
xi_min <- 0.1
xi_max <- -0.8*khi_min**2+0.9
x <- seq(khi_min, khi_max, 0.001)

correl <- data.frame(x=x, y=-0.8*x**2+0.9)
# Create boundary (faster as a dataframe instead of manual geom_lines)
boundary <- data.frame(x=c(khi_min, khi_max, khi_max,khi_min), 
                       y=c(xi_max, xi_max, xi_min, xi_max))

#pdf("mixing.pdf", width=10,height=7)
# A bit of transparency helps
p <- ggplot(data, aes(x=khi, y=xi)) + geom_point(col="red", alpha=1/10)
p <- p + facet_grid(model~type)#, ncol=2)  

p <- p + geom_line(data=correl, aes(x=x, y=y))
p <- p + geom_polygon(data=boundary, aes(x=x, y=y), col="black", alpha=0)
p <- p + coord_cartesian(xlim = c(0, 1.05), ylim = c(0,1)) 
p <- p + xlab(expression(chi)) + ylab(expression(xi))
p <- p + ggtitle("Mixing diagnostic")

p <- p + scale_y_continuous(breaks=seq(0.2,0.8,0.2))
p <- p + scale_x_continuous(breaks=seq(0,1,0.2))
ggsave("mixing_all2.png", width=8, height=10, dpi=500)
