# Impact of slope limitation

library(ggplot2)

folder_pango <- "/wkdir/pae2/praga/parallel/filaments/"
folder <- "/wkdir/pae2/praga/lauritzen/lf/data/"

xlim <- c(0.1,1)
ylim <- c(0,140)
cols <- c("tomato1", "green4")

d <- c()
names <- c( 
    "lf-SLFV-ML-dx0.75-sp-CFL0.8.dat",
    "lf-SBC-dx0.75-sp-CFL1.0.dat",
    "lf-ICON-dx0.77-sp-CFL0.6.dat",
    "lf-pangolin-320lat-sp-CFL0.7.dat")
labels <-  c( 
    "SLFV ML 0.75° sp CFL0.8",
    "SBC 0.75° sp CFL1.0",
    "ICON 0.77° sp CFL0.6",
    "pangolin 0.187° sp CFL0.7")



tau <- c()
lf <- c()
model <- c()
plot_id <- c()
line <- c()

# Append to global vector tau, lf, model, id, line
append_data <- function(name, lmodel, lid, line_style) {
  tmp <- read.table(name)
  tau <<- c(tau, tmp$V1)
  lf <<- c(lf, tmp$V2)
  n <- length(tmp$V1)
  model <<- c(model, rep(lmodel, n))
  plot_id <<- c(plot_id, rep(lid, n))
  line <<- c(line, rep(line_style, n))
}

# Append sp and unlimited data
append_both <- function(name, lmodel, lid) {
  append_data(name, lmodel, id, 1)
  name <- gsub("-sp-", "-un-", name)
  append_data(name, lmodel, id, 2)
}


pdf("limiter_impact_filaments.pdf", width=9, height=4)
# Append data into two array 
# Format : tau, lf, model name, plot id, line style, color
# This will allow for a global (and easier plot)
for (i in seq(length(names)-1)) {
  id <- 1
  if (i > 2) { id <- 2}
  cur <- paste(folder, names[i], sep="")
  append_both(cur, i, id)
}
n <- length(names)
cur <- paste(folder_pango, names[n], sep="")
append_both(cur, n, id)


# Put it into a frame and plot
data <- data.frame(tau, lf, model, plot_id, line)
# Don't forget to factorize
data$model <- factor(data$model)
data$line <- factor(data$line, labels = c("sp", "un"))

# Automatic coloring and line model
p <- ggplot(data,aes(x=tau,y=lf,colour=model)) +
  geom_line(aes(linetype=line), size=0.7)
#data$line <- factor(data$line)
p <- p + facet_wrap(~plot_id, ncol=2)  
p <- p + geom_hline(yintercept=100)

# Customize
p <- p + coord_cartesian(ylim = c(0,140)) 
p <- p + xlab(expression(tau)) + ylab(expression(l[f]))
p <- p + theme(axis.title = element_text(size=20))
## Legend
p <- p + scale_linetype_discrete(name  ="Slope limiter",
                          labels = c("with", "without"))
p <- p + scale_colour_discrete(name  ="Models",
                          labels = labels)
p <- p + ggtitle("Slope limitation impact")
print(p)
dev.off()
