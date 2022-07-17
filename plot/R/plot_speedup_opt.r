# Plot speedup for total advection time with different resolutions only for an
# optimal number of cores (6n^2)


library(ggplot2)
library(reshape2)

# Set the max speedup for all files in directory dir
max_speedup <- function(dir, nb_procs) {
  
  cur <- paste(dir, "/profiling_", ref, sep="")
  time_ref <- read.table(cur)
 
  col <- 6
  t_ref <- max(time_ref[[col]])
 
  # Also reference time with chemistry
  cur<- paste(dir, "/nb_cells", ref, sep="")
  nb_cells <- read.table(cur)[10]
  t_ref_chem <- max(time_ref[[col]] + nb_cells*t_chem_ref)
   
  for (n in nb_procs) {
      cur <- paste(dir, "/profiling_", n, sep="")
      time <- read.table(cur)
      # Only takes total advection time (6th column) and the slowest
      s_loc <- t_ref/max(time[[col]])

      # Add chemistry
      cur<- paste(dir, "/nb_cells", n, sep="")
      # Get total number of interior cells
      nb_cells <- read.table(cur)[10]
      t_chem <- nb_cells*t_chem_ref
      s_loc_chem <- t_ref_chem/(max(time[[col]] + t_chem))
 
      # Update global variables
      smax <<- c(smax, s_loc)
      smax_chem <<- c(smax_chem, s_loc_chem)
 
  }
}

n_max <- 7
# Reference point, either 1 or 3
ref <- 3

# Chemistry, as implemented by Philippe Moinat
t_chem_ref <- 3e-3

nb_procs <- c(1, 6*seq(n_max)**2)
if (ref != 1) {
  nb_procs <- c(3, 6*seq(n_max)**2)
}

smax <- c()
smax_chem <- c()
dir <-"/wkdir/pae2/praga/parallel/speedup/gaussian_1tracer_CFL0.95/320lat_opt"
max_speedup(dir, nb_procs)


# Put it into a frame and plot
data <- data.frame(nb_procs, smax, smax_chem)
# Merge to a single column
data <- melt(data, id="nb_procs")


# Automatic coloring and line model
p <- ggplot(data,aes(x=nb_procs, y=value, colour=variable)) +
geom_point(shape="+", size=15)
p <- p + geom_abline(slope=1./ref, size=1.5) # Reference
p <- p + coord_cartesian(ylim = c(0, 6*n_max**2/ref)) 

#p <- p + theme(axis.title = element_text(size=20))
# Perfect case
p <- p + scale_x_continuous(breaks=c(ref, 6*seq(2, n_max)**2))
p <- p + scale_colour_discrete(name="Type",
    labels=c("Advection", "Advection + chemistry"))

#p <- p + guides(colour=FALSE) # no legend
p <- p + ylab("speedup") + xlab("nb cores")
p <- p + ggtitle("Estimation of chemistry impact (0.28x0.187°)")
# Legend in figure, top left
p <- p + theme(legend.justification=c(0,1), legend.position=c(0,1))
p <- p + theme(text = element_text(size=20))
#print(p)
ggsave("speedup_chem.pdf", width=10,height=10)
