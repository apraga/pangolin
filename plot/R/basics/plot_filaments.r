# List of data files, only for a given resolution
# We substitute the other resolutions later

fname <- "filament2.ps"


# First plot : UCISOM and SBC with same pangolin reference
# Second plot : farsight, homme, MPAS with 3 different pangolin resolutions
names <- c( "lf-UCISOM-CS-dx0.75-sp-CFL1.0.dat",
    "lf-farsight-dx0.75-sp-CFL1.0.dat",
    "lf-SBC-dx0.75-un-CFL1.0.dat",
    "lf-HOMME-dx0.75-sp-p3.dat",
    "lf-MPAS-dx0.75-sp-CFL0.8.dat")
names_pango <-  c("lf-pangolin-138lat-sp-CFL0.7.dat",
    "lf-pangolin-120lat-sp-CFL0.7.dat",
    "lf-pangolin-138lat-sp-CFL0.7.dat",
    "lf-pangolin-40lat-sp-CFL0.7.dat",
    "lf-pangolin-120lat-sp-CFL0.7.dat")

titles <- c("UCISOM, FARSIGHT",
    "SBC, HOMME", "MPAS")

cols <- c('magenta1', 'darkgreen', 'blue', 'orange', 'red')
points <- c(5, 4, 1, 2, 6)

# Add folder names
folder <- "/wkdir/pae2/praga/lauritzen/lf/data/"
folder_pango <- "/wkdir/pae2/praga/parallel/filaments/"
for (i in seq(length(names))) {
  names[i] <- paste(folder, names[i], sep="")
  names_pango[i] <- paste(folder_pango, names_pango[i], sep="")
}

pango_ref <- paste(folder_pango, "lf-pangolin-640lat-sp-CFL0.7.dat", sep="")

# Global plot parameters
x_axis <- c(0:100)
lwidth <- c(1, 2) # Line width
xlim <- c(0.1,1)
ylim <- c(0,140)

#-------------------------------------------------------------------------------
# Plot one filament diagnostic
plot_filament <- function(cur, point, color, style, width=2) {
  a <- read.table(cur)
  x <- unlist(a[1])
  y <- unlist(a[2])
    plot(x,y, type="l",xlim=xlim, ylim=ylim, pch=point,
       col=color, lwd=width, lty=style, xlab="", ylab="", xaxt="n",yaxt="n")
}

# Plot several diagnostics
plot_several <- function(start, end, titl) {
#  if (start == 1) { par(mar=c(2,4,4,2))}
#  else { par(mar=c(5,4,4,2))}
  for (i in seq(start, end)) {
    # Model
    plot_filament(names[i], points[i], cols[i], 1)
    par(new=TRUE)
    # Pangolin with same number of cell
    plot_filament(names_pango[i], 3, cols[i], 2)
    par(new=TRUE)
  }
  # High-res pangolin
  plot_filament(pango_ref, 3, "black", 1)

  # Finalize
  grid()
  abline(100,0)
  title(titl, font.main=1, cex.axis=1.5) # not bold

  # Skip some data 
  if (start != 3) { 
    axis(2, at=seq(0,140,20)) 
    title(ylab=expression(l[f]))
  }
    title(xlab=expression(tau)) 
    axis(1, at=seq(0,1,0.2))
  par(new=FALSE)
}


#-------------------------------------------------------------------------------

#postscript("filament2.ps", height=3, width=9, horizontal=FALSE,
           #onefile=FALSE, paper="special")
           #)#, pointsize=18)
postscript(fname, height=9, width=8.5)#, horizontal=FALSE,
widths <- c(0.5, 0.5, 0.5, 0.5)
heights <- c(3, 3, 4, 6)
l <- layout(matrix(c(1,2,3,4), 2,2, byrow=TRUE), widths=widths)
#heights=heights)#, TRUE)
#par(mfrow=c(2,2))

# Plot several resolutions
par(cex.lab=2)
par(mai=c(0.7,0.7,0.2,0.2))
plot_several(1,2, titles[1])

par(mai=c(0.7,0.2,0.2,0.2))
plot_several(3,4, titles[2])

par(mai=c(0.7,0.7,0.2,0.2))
plot_several(5,5, titles[3])

# Legend in a special (empty) plot
par(new=FALSE)
frame()
par(mar=c(0,0,0,0))
labels <- c("UCISOM-CS sp CFL1.0",
            "farsight sp CFL1.0",
            "SBC un CFL1.0",
            "HOMME sp p3",
            "MPAS sp CFL0.8",
            "Pangolin (0.14x0.9°)",
            "Pangolin",
            "Pangolin",
            "Pangolin",
            "Pangolin",
            "Pangolin")
legend("center", labels,
       col=c(cols, 1, cols), lwd=2, lty= c(rep(1, 5), 1, rep(2,5)),
       ncol=2, bty="n",
       #pch=points, 
       y.intersp=1.5, seg.len=3) # Customization
dev.off()
