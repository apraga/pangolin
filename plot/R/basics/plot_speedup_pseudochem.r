# Add a pseudo chemistry to speedup time. It's just a given time*nb cells
# File format for time: zonal time, advection only zonal time, meridional time,
# meridional advection only time, total advection time
# File format for nb cells : 10 columns, the latest being the total

# Set this for plotting efficiency instead of speepup
savefig <- TRUE#FALSE
titlefig <- "chemistry_speedup.ps"

dir <-"/wkdir/pae2/praga/speedup/90lat_1tracer"
dir_nb <-"/wkdir/pae2/praga/nb_cells/90lat"
n_max <- 63
nb_procs <- seq(3,n_max,3)

# 60 levels, not optimized
#chem_time <- 0.0125
# Philippe version
chem_time <- 3e-3
chem_time_optim <- 0.0125/5

min <- seq(nb_procs)
min_chem <- seq(nb_procs)
min_chem_optim <- seq(length(nb_procs))


cur <- paste(dir, "/profiling_1", sep="")
time_seq <- read.table(cur)
# 90 latitude on an hemisphere
time_seq_chem <- 6*90**2*chem_time
time_seq_chem_optim <- 6*90**2*chem_time_optim

i <- 1
for (n in nb_procs) {
    cur <- paste(dir, "/profiling_", n, sep="")
    time <- read.table(cur)

    cur<- paste(dir_nb, "/nb_cells", n, sep="")
    # Get total number of interior cells
    nb_cells <- read.table(cur)[10]
    time_chem <- nb_cells*chem_time
    time_chem_optim <- nb_cells*chem_time_optim

    time2 <- time
    # Add chemistry time as an offset
    tmp_seq <- time_seq[[6]]
    # Advection time
    min[i] <- min(tmp_seq / time[6])

    j <- 6
    # With pseud chemistry 
    time2[j] <- time[j] + time_chem
    min_chem[i] <- min((tmp_seq + time_seq_chem) / (time2[6]))

    # With pseud chemistry optimized
    time2[j] <- time[j] + time_chem_optim
    min_chem_optim[i] <- min((tmp_seq + time_seq_chem_optim) / (time2[6]))
    i <- i+1
}

extrem <- c(0, n_max)
y <- nb_procs

# Plot min only as the max can be outside the theorical limit
plot_all <- function() {
  # Min only
  k <- 6
  titles <- c("Total zonal", "advec zonal", "total merid", "advec merid",
    "total")
  #X11()

  color <- c("seagreen3", "blue", "red")
  # Advection
  args <- list(nb_procs, min, ylim=extrem, ylab="",xlab="", type="l",
      lty=1)
  args$col <- color[1]
  do.call("plot", args)

  # With chemistry
  par(new=TRUE)
  args[[2]] <- min_chem
  args$col <- color[2]
  do.call("plot", args)

  # With chemistry optim
  par(new=TRUE)
  args[[2]] <- min_chem_optim
  args$col <- color[3]
  do.call("plot", args)

  # Add one tick for first
  axis(1, at=c(3),labels="3")

  title("Estimation of chemistry impact on parallel performance (1x0.67°)")
  legend_pos = "topleft"
  legend(legend_pos, c("advection only", "chemistry(MOCAGE, SUGAR, 60 levels,
  optim)",
  "chemistry (RACMOBUS, Philippe)"), col=color, lty=1)
}


if (savefig) postscript(titlefig,width=10,height=10,pointsize=15)
plot(nb_procs, y, type='l', ylim=extrem, ylab="speedup min", xlab="nb cpus")
par(new=TRUE)
plot_all()
if (savefig) dev.off()
