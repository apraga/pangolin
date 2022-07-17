# Compare filaments diagnostics for models vs pangolin at same number of cells

# We substitute the other resolutions later
library(ggplot2)
first_plot <- FALSE
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
names_pango <-  c("lf-pangolin-70lat-un-CFL1.0.dat",
    "lf-pangolin-120lat-un-CFL1.0.dat",
    "lf-pangolin-124lat-un-CFL1.1.dat",
    "lf-pangolin-139lat-un-CFL1.0.dat", 
    "lf-pangolin-139lat-un-CFL1.0.dat" 
    )


labels <- c( "CLAW un", "FARSIGHT un", "SLFV-ML un", "UCISOM sp", "CAMFV sp", "Pangolin CFL1.0")

tau <- c()
lf <- c()
model <- c()
plot_id <- c()
line <- c()
color <- c()
k <- 1

# Returns the title of the plot (with several models)
get_plot_id <- function(k) {
  if (k <= 2) { id <- "CLAW, FARSIGHT"}
  else if (k <= 4) { id <- "SLFV-ML, UCISOM"}
  else { id <- "CAM-FV"}
  return(id)
}

# Append to global vector tau, lf, model, line, color 
# Plots are grouped in a plot (by plot id)  and by color
# Finally, line style define the line and model is used in the legend
# Color define a "model", that is at a given nb of cells
append_data <- function(name, lcolor, lmodel, lid, line_style) {
  tmp <- read.table(name)
  tau <<- c(tau, tmp$V1)
  lf <<- c(lf, tmp$V2)
  n <- length(tmp$V1)
  model <<- c(model, rep(lmodel, n))
  # Solid line
  line <<- c(line, rep(line_style, n))
  color <<- c(color, rep(lcolor, n))
  plot_id <<- c(plot_id, rep(lid, n))
}

#pdf("filament.pdf", width=12, height=7)
# Append data into two array 
# Format : tau, lf, model name, plot id, line style, color
# This will allow for a global (and easier plot)
for (i in seq(length(names))) {

  id <- get_plot_id(k)
  # Model
  cur <- paste(folder, names[i], sep="")
  append_data(cur, k, labels[k], id, 1)
  
  # Pangolin (same color as model)
  cur <- paste(folder_pango, names_pango[i], sep="")
  append_data(cur, k, labels[k], id, 2)

  # Pangolin high res (in black)
  cur <- paste(folder_pango, "lf-pangolin-640lat-sp-CFL1.0.dat", sep="")
  append_data(cur, 6, "Pangolin", id, 1)

 k <- k + 1
}

# Put it into a frame and plot
data <- data.frame(tau, lf, model, plot_id, line, color)
# Don't forget to factorize
data$line <- factor(data$line)
data$color <- factor(data$color)

# Automatic coloring and line model
p <- ggplot(data,aes(x=tau,y=lf,colour=color)) +
geom_line(aes(linetype=line), size=1.1)
p <- p +facet_wrap(~plot_id, ncol=2)  
p <- p + geom_hline(yintercept=100)

# Customize
p <- p + coord_cartesian(ylim = c(0,140)) 
#p <- p + scale_y_continuous(limits = c(0,140)) 
p <- p + xlab(expression(tau)) + ylab(expression(l[f]))
p <- p + theme(axis.title = element_text(size=20))
# Legend
p <- p + scale_linetype_discrete(name  ="Pangolin",
                          labels = c("0.28x0.18 un°", "same nb of cells un"))
# Palette
cbPalette <- c("darkorange2", "steelblue2", "springgreen4", "orchid3",
"firebrick2", "grey43")

# Legends and color
p <- p + scale_colour_manual( name  ="Models", labels = labels, 
    values=cbPalette)

p <- p + ggtitle("Filaments diagnostic comparison")
print(p)
#dev.off()
