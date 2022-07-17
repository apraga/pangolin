# Plot speedup for total advection time with different resolutions for every
# number of cores


library(ggplot2)
library(reshape2)

# Returns the max speedup for all files in directory dir
max_speedup <- function(dir, nb_procs) {
  
  smax <- c()
  
  cur <- paste(dir, "/profiling_", ref, sep="")
  time_ref <- read.table(cur)
  
  col <- 6
  t_ref <- max(time_ref[[col]])

  for (n in nb_procs) {
      cur <- paste(dir, "/profiling_", n, sep="")
      time <- read.table(cur)
      # Only takes total advection time (6th column) and the slowest
      smax <- c(smax, t_ref/max(time[[col]]))
  }
  return(smax)
}

n_max <- 126
# Reference point, either 1 or 3
ref <- 3

nb_procs <- c(1, seq(3,n_max,3))
if (ref != 1) {
  nb_procs <- seq(3, n_max,3)
}

dir <-"/wkdir/pae2/praga/parallel/speedup/gaussian_1tracer_CFL0.95/80lat_small"
s_80 <- max_speedup(dir, nb_procs)
dir2 <- gsub("80", "160", dir)
s_160 <- max_speedup(dir2, nb_procs)
dir2 <- gsub("80", "320", dir)
s_320 <- max_speedup(dir2, nb_procs)

# Put it into a frame and plot
data <- data.frame(nb_procs, s_80, s_160, s_320)
# Merge to a single column
data <- melt(data, id="nb_procs")


# Automatic coloring and line model
p <- ggplot(data,aes(x=nb_procs, y=value, colour=variable)) + 
geom_line(size=1) + geom_point(size=3)
p <- p + geom_abline(slope=1./ref, size=1.5) # Reference
p <- p + coord_cartesian(ylim = c(0, n_max/ref)) 

#p <- p + theme(axis.title = element_text(size=20))
# Perfect case
p <- p + scale_x_continuous(breaks=c(ref, 6*seq(4)**2))
p <- p + scale_colour_discrete(name="Nb lat/2",
    labels=c("1.13x0.75°","0.56x0.38°","0.28x0.19°"))

p <- p + ylab("speedup") + xlab("nb cores")
p <- p + ggtitle("Resolution impact on speedup (1 tracer)")
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1),
               text = element_text(size=20))
#print(p)
ggsave("speedup_resolution.pdf", width=10,height=10)
