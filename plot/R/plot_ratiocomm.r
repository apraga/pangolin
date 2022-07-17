# Plot max, min, and mean of different advections times.
# File format : proc zonal_total zonal_advec merid_total merid_advec advection 

dir <- "/wkdir/pae2/praga/speedup/90lat"

nb_procs <- seq(3,63,3)

# We store only ratio communications/advection
nb_cols = 2
#mean <- matrix(nrow=length(nb_procs), ncol=nb_cols)
#min <- matrix(nrow=length(nb_procs), ncol=nb_cols)
ratio <- matrix(nrow=length(nb_procs), ncol=nb_cols)

i <- 1
for (n in nb_procs) {
    cur <- paste(dir, "/profiling_", n, sep="")
    time <- read.table(cur)
    # Communication times : zonal, then meridional
    index <- which.max(time[,2])
    ratio[i,1] <- (time[index, 2]-time[index, 3])/time[index, 2]

    index <- which.max(time[,4])
    ratio[i,2] <- (time[index, 4]-time[index, 5])/time[index, 4]
    i <- i+1
}
#

#postscript("ratio_comm.ps",width=10,height=10,pointsize=15)
## Space for x title
##par(mfrow=c(1,2))
extrem <- c(min(ratio[,]), max(ratio[,]))
plot(nb_procs, ratio[,1], type='l', ylim=extrem, lty=1, ylab="",col="red")
par(new=TRUE)
plot(nb_procs, ratio[,2], type='l', ylim=extrem, lty=2, ylab="",col="blue")
#plot(nb_procs, max[,1], type='l', ylim=extrem, lty=2, ylab="",col="red")
#par(new=TRUE)
#plot(nb_procs, mean[,1], type='l', ylim=extrem, lty=1, ylab="",col="red")
#
#par(new=TRUE)
#plot(nb_procs, max[,2], type='l', ylim=extrem, lty=2, ylab="",col="blue")
#par(new=TRUE)
#plot(nb_procs, mean[,2], type='l', ylim=extrem, lty=1, ylab="",col="blue")
#mtext("zonal", side=2, padj=-4, col="red", cex=1.2)
#mtext("meridional", side=4, padj=1, col="blue", cex=1.2)
#
#title("Communication impact (1x0.67° at the equator)")
#legend("topright", c("max", "mean"), lty=c(2,1))
#dev.off()
