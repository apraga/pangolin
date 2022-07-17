# List of data files, only for a given resolution
# We substitute the other resolutions later

fname <- "filament2.ps"

# First plot
names <- c( "lf-UCISOM-CS-dx0.75-sp-CFL1.0.dat",
    "lf-SBC-dx0.75-un-CFL1.0.dat",
    "lf-farsight-dx0.75-sp-CFL1.0.dat",
    "lf-HOMME-dx0.75-sp-p3.dat",
    "lf-MPAS-dx0.75-sp-CFL0.8.dat"
    );

folder <- "/wkdir/pae2/praga/lauritzen/lf/data/"
for (i in seq(length(names))) {
  names[i] <- paste(folder, names[i], sep="")
}

folder2 <- "/wkdir/pae2/praga/parallel/filaments/"
names_pango <- c()
i <- 1
for (n in c('40','80', '120')) {
  names_pango[i] <-  paste("lf-pangolin-",n, "lat-sp-CFL0.7.dat", sep='')
  names_pango[i] <- paste(folder2, names_pango[i], sep="")
  i <- i+1
}
cols <- c('orange', 'green2', 'blue', 'purple', 'grey', 'black')
points <- c(4, 15, 22)

# Colors and point styles

# Global plot parameters
x_axis <- c(0:100)
lwidth <- c(1, 2) # Line width
xlim <- c(0.1,1)
ylim <- c(0,140)

#-------------------------------------------------------------------------------
# Plot one filament diagnostic
plot_filament <- function(cur, point, color, linew, style) {
  a <- read.table(cur)
  x <- unlist(a[1])
  y <- unlist(a[2])
    plot(x,y, type="o",xlim=xlim, ylim=ylim, col=color, pch=point,
    lty=style, lwd=linew, xlab="", ylab="", xaxt="n",yaxt="n")
}

# Plot all diagnostic for 1 resolution
for (i in seq(length(names))) {
  plot_filament(names[i], points[i],cols[i], 2, 1);
  par(new=TRUE)
  #plot_filament(names_pango[i], points[i],cols[i], 3, 2);
  #par(new=TRUE)
}
axis(1, at=seq(0,1,0.2))
axis(2, at=seq(0,140,20))
#title(paste(res, "unlimited and shape-preserving"),
title("lol",
    font.main=1, # not bold
    xlab=expression(tau), ylab=expression(l[f])) #subscript
