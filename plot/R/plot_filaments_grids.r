# Compare filaments diagnostics for models vs pangolin high-res

library(ggplot2)
fname <- "filament2.ps"

folder_pango <- "/wkdir/pae2/praga/parallel/filaments/"
folder <- "/wkdir/pae2/praga/lauritzen/lf/data/"

xlim <- c(0.1,1)
ylim <- c(0,140)
cols <- c("tomato1", "green4")

d <- c()
names <- c( "lf-CLAW-dx0.75-un-CFL0.95.dat",
    "lf-farsight-dx0.75-un-CFL1.0.dat",
    "lf-SLFV-ML-dx0.75-un-CFL0.8.dat",
    "lf-UCISOM-dx0.75-sp-CFL5.5.dat", 
    "lf-camfv-dx0.75-sp-CFL1.2.dat" 
    )

names_pango <-  c("lf-pangolin-139lat-un-CFL1.0.dat",
    "lf-pangolin-640lat-sp-CFL1.0.dat")


labels <- c( "CLAW CFL0.95 un", "FARSIGHT CFL1.0", "SLFV-ML un CFL0.8", "UCISOM CFL1.0", "CAMFV CFL1.2", "Pangolin 0.65x0.43° CFL1.0", "Pangolin 0.14x0.09° CFL1.0")

tau <- c()
lf <- c()
plot_id <- c()
line <- c()
color <- c()
k <- 1

# Returns the title of the plot (with several models)
get_plot_id <- function(k) {
  if (k < 4) { id <- "Non-lat-lon grid"}
  else { id <- "Lat-lon grid"}
  return(id)
}

# Append to global vector tau, lf, model, line, color 
# Plots are grouped in a plot (by plot id)  and by color
# Finally, line style define the line and model is used in the legend
# Color define a "model", that is at a given nb of cells
append_data <- function(name, lcolor, lid, line_style) {
  tmp <- read.table(name)
  tau <<- c(tau, tmp$V1)
  lf <<- c(lf, tmp$V2)
  n <- length(tmp$V1)
  # Solid line
  line <<- c(line, rep(line_style, n))
  color <<- c(color, rep(lcolor, n))
  plot_id <<- c(plot_id, rep(lid, n))
}

pdf("filaments_compare.pdf", width=10, height=4)
# Append data into two array 
# Format : tau, lf,  plot id, line style, color
# This will allow for a global (and easier plot)
for (i in seq(length(names))) {

  id <- get_plot_id(k)
  # Model
  cur <- paste(folder, names[i], sep="")
  append_data(cur, k, id, 2)

  k <- k + 1
}

# Pangolin data
cur <- paste(folder_pango, names_pango[1], sep="")
id <- get_plot_id(1)
append_data(cur, k, id, 1)
  k <- k + 1
cur <- paste(folder_pango, names_pango[2], sep="")
id <- get_plot_id(4)
append_data(cur, k,  id, 1)

# Put it into a frame and plot
data <- data.frame(tau, lf, plot_id, line, color)
# Don't forget to factorize
data$line <- factor(data$line)
data$color <- factor(data$color)


# Automatic coloring and line model
p <- ggplot(data,aes(x=tau,y=lf)) +
  geom_line(aes(linetype=line, colour=color))#, size=1)

p <- p + facet_wrap(~plot_id, ncol=2)  
# Custom facet labels
p <- p + theme(strip.text.x = element_text(size=12, face="bold"))


p <- p + geom_hline(yintercept=100)

# Customize
p <- p + coord_cartesian(ylim = c(0,140)) 
p <- p + xlab(expression(tau)) + ylab(expression(l[f]))
p <- p + theme(axis.title = element_text(size=20))
# Palette
cbPalette <- c("darkorange2", "steelblue2", "springgreen4", "orchid3",
"firebrick2", "grey43", "black")

# Legends and color
p <- p + scale_colour_manual( name  ="Models", labels = labels, 
   values=cbPalette)
p <- p + scale_linetype_discrete( guide=FALSE)

print(p)
dev.off()
