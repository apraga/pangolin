# Plot speedup/efficiency from different advection times (max, min, mean)
# File format : zonal time, advection only zonal time, meridional time,
# meridional advection only time, total advection time

# Set this for plotting efficiency instead of speepup
is_efficiency <- FALSE
savefig <- TRUE#FALSE

dir <-"/wkdir/pae2/praga/speedup/90lat"
n_max = 63
nb_procs <- seq(3,n_max,3)

mean <- matrix(nrow=length(nb_procs), ncol=5)
min <- matrix(nrow=length(nb_procs), ncol=5)
max <- matrix(nrow=length(nb_procs), ncol=5)

cur <- paste(dir, "/profiling_1", sep="")
time_seq <- read.table(cur)

i <- 1
for (n in nb_procs) {
    cur <- paste(dir, "/profiling_", n, sep="")
    time <- read.table(cur)
    for (j in seq(2,6)) {
      time[j] <- time_seq[[j]] / time[j] 
      mean[i,j-1] <- colMeans(time[j])
      max[i,j-1] <- max(time[j])
      min[i,j-1] <- min(time[j])
    }
    i <- i+1
}

extrem <- c(0, n_max)
y <- nb_procs

if (is_efficiency) {
  extrem <- c(0,1)
  mean <- mean/nb_procs
  max <- max/nb_procs
  min <- min/nb_procs
  y <- rep(1, length(nb_procs))
}

# Plot max, mean, min
plot_all <- function(k) {
  titles <- c("Total zonal", "advec zonal", "total merid", "advec merid",
      "total")
    #X11()
    args <- list(nb_procs, mean[,k], ylim=extrem, ylab="",xlab="", type="o",
        lty=1)
    do.call("plot", args)

    par(new=TRUE)
    args[[2]] <- max[,k]
    args$lty=2
    do.call("plot", args)

    par(new=TRUE)
    args[[2]] <- min[,k]
    args$lty=3
    do.call("plot", args)

    # Add one tick for first
    axis(1, at=c(3),labels="3")


    title("Parallel performance (1x0.67° at the Equator)")
    legend_pos = "topleft"
    if (is_efficiency) legend_pos = "bottomright"
    legend(legend_pos, c("mean", "max", "min"), lty=c(1,2,3))
}


if (savefig) postscript("speedup.ps",width=10,height=10,pointsize=15)
plot(nb_procs, y, type='l', ylim=extrem,ylab="speedup", xlab="nb cpus")
par(new=TRUE)
plot_all(5)
if (savefig) dev.off()
