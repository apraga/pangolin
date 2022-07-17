# Compare filaments diagnostics for models vs pangolin high-res

library(ggplot2)

folder <- "/wkdir/pae2/praga/parallel/filaments/"

xlim <- c(0.1,1)
ylim <- c(0,140)
cols <- c("tomato1", "green4")

d <- c()
names <-  c("lf-pangolin-320lat-sp-CFL0.7.dat",
    "lf-pangolin-320lat-sp-nocorr-CFL0.7.dat")

tau <- c()
lf <- c()
color <- c()
k <- 1

# Append to global vector tau, lf, model, line, color 
append_data <- function(name, lcolor) {
  tmp <- read.table(name)
  tau <<- c(tau, tmp$V1)
  lf <<- c(lf, tmp$V2)
  n <- length(tmp$V1)
  color <<- c(color, rep(lcolor, n))
}

# Note : if "q" appears in pdf, see ?pdf help 
pdf("windscorr_impact_filaments.pdf", width=9, height=5)
for (cur in names) {
  cur <- paste(folder, cur, sep="")
  append_data(cur, k)
 k <- k + 1
}

# Put it into a frame and plot
data <- data.frame(tau, lf, color)
# Factor for discrete scale and custom legend
data$color <- factor(data$color)


# Automatic coloring and line model
p <- ggplot(data,aes(x=tau,y=lf,colour=color)) +
  geom_line(size=1) + geom_point(aes(shape=color), size=5)

p <- p + geom_hline(yintercept=100)

# Customize
p <- p + scale_shape_manual(name="Winds correction", labels=c("with",
"without"), values=c(25,20))

p <- p + coord_cartesian(ylim = c(0,140)) 
p <- p + xlab(expression(tau)) + ylab(expression(l[f]))
p <- p + theme(axis.title = element_text(size=20))
p <- p + ggtitle("Winds correction impact (sp, 0.28x0.18°)")
p <- p + guides(colour=FALSE)

print(p)
dev.off()
