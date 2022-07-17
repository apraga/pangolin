# Compare cv rates vs nb cells for different models
# This version plots for limited scheme, with hourdin comparison

library(ggplot2)
library(scales) # for log

#-------------------------------------------------------------------------------
# Nb cells from resolution
latlon <- function(lon) { return(floor(360/lon)*floor(180/lon)) }

homme <- function(lon) { return(6*(90/(3*lon))**2) }

farsight <- function(lon) { return(6*(90/lon)**2) }

pango <- function(lon) { 
  res <- floor(120/lon)
  res <- floor(0.5*(res+1))
  res <- 6*res^2
  return(res) 
}

mpas <- function(lon) { 
  res <- pango(lon)
  return(2.24*res) 
}


# Append to global vectors
append_data <- function(dat, is_pango, nbcells, nmodel, ntype) {
  nb_cells <<- c(nb_cells, nbcells, nbcells)
  size <- length(dat$V2)
  # Pangolin has one more column
  if (is_pango) { 
    err <<- c(err, dat$V3, dat$V4)
  } else {
    err <<- c(err, dat$V2, dat$V3)
  }
  err_type <<- c(err_type, rep("l2", size), rep("loo", size))
  model <<- c(model, rep(nmodel, 2*size))
  type <<- c(type, rep(ntype, 2*size))
}

#-------------------------------------------------------------------------------

pdf("cv_rate_compare_hourdin.pdf", width=10,height=7)
folder_pango <- "/wkdir/pae2/praga/parallel/cv_rate/"
folder <- "/wkdir/pae2/praga/lauritzen/conv/data/"
labels <- c("UCISOM CFL1.0", "FARSIGHT CFL1.0", "SBC CFL1.0", "HOMME p3",
            "MPAS CFL0.8")


names <- c( "conv-UCISOM-sp-gs-CFL0.8.dat",
           "conv-FARSIGHT-sp-gs-CFL1.0.dat",
           "conv-SBC-sp-gs-CFL1.0.dat",
           "conv-HOMME-sp-gs-p3.dat",
           "conv-MPAS-sp-gs-CFL0.8.dat")


nb_cells <- c()
err <- c()
err_type <- c()
model <- c()
type <- c()

# Plot with Hourdin test, for limited versio
##both unlimited and limited
for (ttype in c("sp", "hourdin-rot", "hourdin-norot")){#, "un")) {
  name <- "conv-pangolin-sp-CFL0.7.dat"
  name <- gsub("sp", ttype, name)

  dat_pango <- read.table(paste(folder_pango, name, sep=""))
  append_data(dat_pango, TRUE, dat_pango$V2, "PANGOLIN CFL0.7", ttype)

  for (i in seq(length(names))) {
    #cur <- gsub("sp", ttype, names[i])
    cur <- names[i]
    cur <- paste(folder, cur, sep="")
    # Skip missing data
    if (!file.exists(cur)) { 
    print(cur)
    next }
    dat <- read.table(cur)

    if (i == 1 || i == 3) { n <- latlon(dat$V1) }
    else if ( i == 2) { n <- farsight(dat$V1) }
    else if ( i == 4) { n <- homme(dat$V1) }
    else { n <- mpas(dat$V1)} # estimation

    append_data(dat, FALSE, n, labels[i], ttype)
  }
}
data <- data.frame(nb_cells, err, err_type, model, type)
p <- ggplot(data, aes(x=nb_cells, y=err, colour=model)) + geom_line(aes(linetype=type))
p <- p + facet_wrap(~err_type, ncol=2)  
p <- p + xlab("nb cells") + ylab("error")
#p <- p + theme(axis.title = element_text(size=18))
p <- p + ggtitle("Sp convergence rates")

# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(1e-1,1e-2,1e-3,1e-4,1e-5))
p <- p + scale_x_continuous(trans=log2_trans(),
                            breaks=c(1e3, 1e4, 1e5,1e6))
print(p)
dev.off()
