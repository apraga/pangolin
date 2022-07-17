# Plot several speedup/efficiency from different advection times (mean only)
# File format : zonal time, advection only zonal time, meridional time,
# meridional advection only time, total advection time
# Arguments : directory

is_efficiency <- TRUE#FALSE
savefig <- TRUE#FALSE


plot_speedup <- function(nb_procs, dir) {
  
  
  mean <- matrix(nrow=length(nb_procs), ncol=5)
  
  cur <- paste(dir, "/profiling_1", sep="")
  time_seq <- read.table(cur)
  
  i <- 1
  for (n in nb_procs) {
      cur <- paste(dir, "/profiling_", n, sep="")
      time <- read.table(cur)
      for (j in seq(2,6)) {
        time[j] <- time_seq[[j]] / time[j] 
        mean[i,j-1] <- colMeans(time[j])
      }
      i <- i+1
  }

  if(is_efficiency) mean <- mean/nb_procs
  return(mean)
}

n_max = 63
nb_procs <- seq(3,n_max,3)

mean90 <- plot_speedup(nb_procs, "/wkdir/pae2/praga/speedup/90lat")
mean180 <- plot_speedup(nb_procs, "/wkdir/pae2/praga/speedup/180lat")

y <- nb_procs
if (is_efficiency) y <- rep(1, length(nb_procs))
extrem <- c(0, max(y))
args <- list(nb_procs, y, xlab="", ylab="", type="l",ylim=extrem,col="black")

if (savefig) postscript("bothspeedup.ps",width=3,height=9,pointsize=15)
do.call("plot", args)

par(new=TRUE)
args[[2]] <- mean90[,5]
args$col <- "blue"
args$type <- "o"
do.call("plot", args)

par(new=TRUE)
args[[2]] <- mean180[,5]
args$col <- "red"
args$ylab <- "mean(speedup)"
if (is_efficiency) args$ylab <- "mean(efficiency)"
args$xlab <- "nb cpus"
do.call("plot", args)

title_s <- "Speedup vs resolution at the Equator"
if (is_efficiency) title_s <- "Efficiency vs resolution at the Equator"
title(title_s)

legend_pos = "topleft"
if (is_efficiency) legend_pos = "bottomright"
legend(legend_pos, c("1x0.67°", "0.5x0.33°"), col=c("blue", "red"), lwd=1)

if (savefig) dev.off()
