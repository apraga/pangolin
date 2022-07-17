# First columun is l_2, secnode is l_oo
folder <- "/wkdir/pae2/praga/lauritzen/conv/data/"

first_plot <- FALSE
fname <- "cv_rate2.ps"

#if (first_plot) {
if (first_plot) {
  names <- c("conv-UCISOM-sp-gs-CFL0.8.dat",
             "conv-CAM-FV-sp-gs-CFL0.2.dat",
             "conv-pangolin-sp-CFL0.7.dat")
  cols <- c('magenta1', 'green2', 'blue')
} else {
  names <- c("conv-CLAW-sp-gs-CFL0.95.dat",
            "conv-FARSIGHT-sp-gs-CFL1.0.dat",
            "conv-SLFV-ML-sp-gs-CFL0.8.dat")
cols <- c('orange', 'red', 'blue')
}

points <- c(1, 15, 2, 17)

# Revert x order
xlim <- c(3.0, 0.1875)
ylim <- c(1e-5,10)

xlabel <- c(3, 1.5, 0.75, 0.375, 0.1895 )
ylabel <- c('10^-5', '10^-4', '10^-3', '10^-2', '10^-1', '10^0' ,
      '10^1')
ylabel <- parse(text=ylabel)

#-------------------------------------------------------------------------------
# Norm is 2 for 2-norm, otherwise third)
        plot_error <- function(name, col, pch, norm) {
  c <- read.table(name)
  x <- unlist(c[1])
  if (norm == 2) { y <- unlist(c[2]) } 
  else { y <- unlist(c[3])}
  plot(x,y, log="xy", type="o",xlim=xlim, ylim=ylim, 
      xaxt="n", yaxt="n", xlab="", ylab="", col=col, pch=pch)

}

# Plot errors for a given norm and case (unlimited/limited)
# case= un, sp
# norm = 2, oo
plot_all_error <- function(norm, case) {
  par(mar=c(2,4,2,0))
  for (i in seq(1, length(names))) {
    cur <- paste(folder, names[i], sep="")
    cur <- gsub('sp', case, cur)
    if (file.exists(cur)) {
      plot_error(cur, cols[i],  points[i], 2)
      par(new=TRUE)
    }
  }

  # Plot first and second order
  args <- list(xlabel, 0.5*xlabel, type="l", log="xy", xlim=xlim, ylim=ylim, 
      xaxt="n", yaxt="n", xlab="", ylab="")
  do.call("plot", args)
  par(new=TRUE)
  args[[2]] <- 0.0015*xlabel^2
  do.call("plot", args)

  # Custom axis and grid
  axis(1, at=xlabel, labels=xlabel)
  axis(2, at=10^(-5:1), labels=ylabel, las=TRUE)
  abline(h=10^(-5:1), v=xlabel, col="lightgray",lty=3)
  par(new=FALSE)

  # Custom axis title
  if (case == "un") {
    if (norm == 2) { title(ylab=expression(l[2]),cex.lab=1.5) }
    else { title(ylab=expression(l[infinity]),cex.lab=1.5) }
  }
  if (norm == 2) { 
    if (case == "un") { title("Unlimited") }
    else { title("Shape preserving") }
  }
}
#-------------------------------------------------------------------------------


postscript(fname, height=6, width=9)#, horizontal=FALSE,
widths <- c(0.4, 0.4, 0.2)
layout(matrix(c(1,2,3,4, 5, 0), 2, 3, byrow=FALSE), width=widths)


plot_all_error(2, 'un')
plot_all_error(3, 'un')
plot_all_error(2, 'sp')
plot_all_error(3, 'sp')

# Legend in a special (empty) plot
par(new=FALSE)
frame()
par(mar=c(0,0,0,0))
if (first_plot) {
  labels <- c("UCISOM CN1.0", "CAMFV CFL1.2", "PANGOLIN") 
} else {
  labels <- c("CLAW CFL0.95", "FARSIGHT CFL1.0", "SLFV-ML CFL0.8") 
}
legend("center", labels,
       col=cols, lwd=1, pch=points, 
       y.intersp=1.5, seg.len=3) # Customization
dev.off()
