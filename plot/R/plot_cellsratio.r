# Plot ratio interior cells/halo
# File format : interior/halo cells zonal flux, interior/halo gradient for
# meridional, interior/hal for meridional flux, interior/halo for meridional
# slope

dir <- "/wkdir/pae2/praga/nb_cells/90lat"

nb_procs <- seq(3,63,3)

nb_cols = 5
max <- matrix(nrow=length(nb_procs), ncol=nb_cols)
min <- matrix(nrow=length(nb_procs), ncol=nb_cols)
mean <- matrix(nrow=length(nb_procs), ncol=nb_cols)

i <- 1
for (n in nb_procs) {
    cur <- paste(dir, "/nb_cells", n, sep="")
    nb <- read.table(cur)
    tmp <- matrix(nrow=dim(nb)[1], ncol=nb_cols)
    for (j in seq(nb_cols)) {
      # Last colomun is the total number of cells
      if (j == nb_cols) tmp[,j] <- nb[,2*j]
      else {
        for (k in seq(1,dim(nb)[1])) {
          tmp[k,j] <- nb[k,2*j]/nb[k,2*j+1]
        }
      }
      max[i,j] <- max(tmp[,j])
      min[i,j] <- min(tmp[,j])
      mean[i,j] <- mean(tmp[,j])
    }
    i <- i+1
}

titles <- c("Flux z", "Grad m", "flux m", "slope m")

postscript("ratio_nbcells.ps",width=10,height=10,pointsize=15)
  extrem <- c(min(max[,1:3],mean[,1:3],min[,1:3]), max(max[,1:3], mean[,1:3],min[,1:3]))
  plot(nb_procs, max[,1], type="o", ylim=extrem, ylab="")
  par(new=TRUE)
  plot(nb_procs, max[,3], col="red", type="o", ylim=extrem,
      ylab="max(#interior/#halo)")
  legend("bottomleft", c("Zonal fluxes", "Meridional fluxes"),
      col=c("black","red"),lwd=2)
  # Add one tick for first
  axis(1, at=c(3),labels="3")
  title("Impact of halo cells")
  dev.off()

## Plot everything
#par(mfrow=c(2,2))
#for (k in seq(1,4)) {
#  extrem <- c(min(max[,k],mean[,k],min[,k]), max(max[,k], mean[,k],min[,k]))
#
#  plot(nb_procs, mean[,k], lty=1, type="o", ylim=extrem, ylab="")
#  par(new=TRUE)
#  plot(nb_procs, max[,k], lty=2, type="o", ylim=extrem, ylab="")
#  par(new=TRUE)
#  plot(nb_procs, min[,k], lty=3, type="o", ylim=extrem, ylab="")
##  par(new=TRUE)
##  plot(nb_procs, log(mean[,5]), type="o", ylim=extrem, ylab="",col="red")
#  title(titles[k])
#  legend("topright", c("mean", "max", "min"), lty=c(1,2,3))
#  }
