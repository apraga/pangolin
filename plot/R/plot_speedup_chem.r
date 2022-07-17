# Plot speedup for total advection time with different chemistry times

library(ggplot2)
library(reshape2)

# Returns the max speedup for all files in directory dir
max_speedup <- function(dir, nb_procs) {
  
  smax <- c()
  
  cur <- paste(dir, "/profiling_", ref, sep="")
  time_ref <- read.table(cur)
  
  for (n in nb_procs) {
      cur <- paste(dir, "/profiling_", n, sep="")
      time <- read.table(cur)
      for (i in seq(2,length(time_ref))) {
        time[[i]] <- time_ref[[i]]/time[[i]]
      }
      # Only takes total advection time (6th column)
      smax <- c(smax, max(time[[6]]))
  }
  return(smax)
}

n_max <- 126
# Reference point, either 1 or 3
ref <- 3

nb_procs <- c(1, seq(3,n_max,3))
if (ref != 1) {
  nb_procs <- seq(3,n_max,3)
}

dir <-"/wkdir/pae2/praga/parallel/speedup/320lat_1tracer_CFL0.7"
s_320 <- max_speedup(dir2, nb_procs)

# Put it into a frame and plot
data <- data.frame(nb_procs, s_80, s_160, s_320)
# Merge to a single column
data <- melt(data, id="nb_procs")

pdf("speedup_chem.pdf", width=10,height=7)

# Automatic coloring and line model
p <- ggplot(data,aes(x=nb_procs, y=value, colour=variable)) + geom_line() + geom_point()
p <- p + geom_abline(slope=1./ref) # Reference
p <- p + coord_cartesian(ylim = c(0, n_max/ref)) 

#p <- p + theme(axis.title = element_text(size=20))
# Perfect case
p <- p + scale_x_continuous(breaks=c(ref, 6*seq(4)**2))
p <- p + scale_colour_discrete(name="Nb lat/2",
    labels=c("80","160","320"))

#p <- p + guides(colour=FALSE) # no legend
p <- p + ylab("mean") + xlab("nb cores")
p <- p + ggtitle("Resolution impact on speedup (1 tracer)")
print(p)
dev.off()
