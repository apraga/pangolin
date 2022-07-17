# Show how the mininal resolution is computed
library(ggplot2)
library(scales) # for log

folder_pango <- "/wkdir/pae2/praga/parallel/min_res/"

# Plot only unlimited
name <- "min-res-pangolin-un-CFL1.0.dat"

dat <- read.table(paste(folder_pango, name, sep=""))
x <- c(dat$V1)
err <- c(dat$V3)

#Linear regression
reg <- lm(log(err)~log2(x))
# Reconstruct line
x2 <- seq(0.3, 0.03, -0.01)
a <- as.numeric(reg$coefficients[2])
b <- as.numeric(reg$coefficients[1])
# Log is applied later
err2 <- exp(a*log2(x2))*exp(b)

intersect <- exp(log(2)*(log(0.033)-b)/a)

pdf("min_res_estim.pdf", width=10, height=8)
data <- data.frame(x, err)
data2 <- data.frame(x2, err2)

p <- ggplot(data, aes(x=x, y=err)) 
p <- p + geom_line(colour="blue") + geom_point(colour="blue")
# Add custom line
p <- p + geom_line(data=data2, aes(x=x2, y=err2), linetype="dashed")
p <- p + xlab(expression(paste(Delta,lambda))) + ylab(expression(l[2]))

p <- p + geom_hline(yintercept = 0.033)
p <- p + ggtitle("Minimal resolution estimation")
#p <- p + theme(axis.title = element_text(size=18))

# Log scale
# Custom reverse log
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
        log_breaks(base = base), 
        domain = c(1e-100, Inf))
}
# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(1e-5, 1e-4, 1e-3, 0.033, 1e-2,1e-1,1))
p <- p + scale_x_continuous(trans=reverselog_trans(2),
                            breaks=c(0.3761, 0.0938, 0.1877,0.0938, intersect))
p <- p + coord_cartesian(ylim = c(1e-3, 1))                            
#p <- p + ggtitle("Pole impact on convergence rate (Hourdin, one period)")
p <- p + theme(axis.title = element_text(size=22))
print(p)
dev.off()
