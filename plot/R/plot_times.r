# Plot max, min, and mean of different advections times.
# File format : proc zonal_total zonal_advec merid_total merid_advec advection 

dir <- "/wkdir/pae2/praga/speedup/90lat"
#dir <- "../../data/speedup/90lat"

nb_procs <- seq(3,63,3)

# We store communications time also. Format :
# total_z advec_z comm_z total_m advec_m comm_m total
nb_cols = 7
mean <- matrix(nrow=length(nb_procs), ncol=nb_cols)
min <- matrix(nrow=length(nb_procs), ncol=nb_cols)
max <- matrix(nrow=length(nb_procs), ncol=nb_cols)

i <- 1

for (n in nb_procs) {
    cur <- paste(dir, "/profiling_", n, sep="")
    time <- read.table(cur)
    # total_z advec_z comm_z
    k <- 2
    for (j in seq(1,nb_cols)) {
      vec <- time[k]
      # Communications are deduced
      if (j == 3 || j == 6){
        k <- k-1
        vec <- time[k-1]-time[k]
      }

      mean[i,j] <- colMeans(vec)
      max[i,j] <- max(vec)
      min[i,j] <- min(vec)
       k <- k + 1
    }

    i <- i+1
}


plot_time <- function(k) {
# Space for x title
  titles <- c("Total zonal", "advec zonal", "comm zonal", "total merid", "advec merid",
      "comm merid")
    #X11()
    extrem = c(min(mean[,k], max[,k], mean[,k]), max(mean[,k], max[,k], mean[,k]))
    args <- list(nb_procs, mean[,k], type='l', ylim=extrem, lty=1, ylab="",xlab="")
    do.call("plot", args)

    par(new=TRUE)
    args[[2]] <- max[,k]
    args$lty = 2
    do.call("plot", args)

    par(new=TRUE)
    args[[2]] <- min[,k]
    args$lty = 3
    do.call("plot", args)

    # Add one tick for first
    axis(1, at=c(3),labels="3")

    title(titles[k])
    legend("topright", c("mean", "max", "min"), lty=c(1,2,3))
}

postscript("alltimes.ps",width=10,height=10,pointsize=15)
# Plot all
par(mfrow=c(2,2), mar=c(2,2,2,2))
#par(mfrow=c(2,3), mar=c(2,2,2,2))
  plot_time(2)
  plot_time(3)
  plot_time(5)
  plot_time(6)

  plot_time(1)
  plot_time(4)
dev.off()
