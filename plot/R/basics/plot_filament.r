# For quick plotting

# List of data files, only for a given resolution
# We substitute the other resolutions later

first_plot <- FALSE
fname <- "filament2.ps"

folder <- "/wkdir/pae2/praga/parallel/filaments/"

#postscript(fname, height=3, width=9)#, horizontal=FALSE,
xlim <- c(0.1,1)
ylim <- c(0,140)

# Plot one filament diagnostic
for (lat in c(640,1280)) {
  cur <- paste(folder, "lf-pangolin-",lat,"lat-sp-CFL0.7.dat", sep="")
  a <- read.table(cur)
  x <- unlist(a[1])
  y <- unlist(a[2])
  plot(x,y, type="o" ,xlim=xlim, ylim=ylim)
  par(new=TRUE)
}

abline(100,0)

#dev.off()
