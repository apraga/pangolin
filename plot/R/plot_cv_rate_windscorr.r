# Slope limitation impact

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
append_data <- function(dat, nbcells, ntype, color) {
  nb_cells <<- c(nb_cells, nbcells, nbcells)
  size <- length(dat$V2)
  # Pangolin has one more column
  err <<- c(err, dat$V3, dat$V4)
  err_type <<- c(err_type, rep("l2", size), rep("loo", size))
  type <<- c(type, rep(ntype, 2*size))
  colors <<- c(colors, rep(color, size))
  colors <<- c(colors, rep(color+1, size))
}

#-------------------------------------------------------------------------------

# Note : if "q" appears in pdf, see ?pdf help 
pdf("windscorr_impact_cv_rate.pdf", width=9,height=5)
folder_pango <- "/wkdir/pae2/praga/parallel/cv_rate/"
folder <- "/wkdir/pae2/praga/lauritzen/conv/data/"
labels <- c("UCISOM CFL1.0", "FARSIGHT CFL1.0", "SBC CFL1.0", "HOMME p3",
            "MPAS CFL0.8")

nb_cells <- c()
err <- c()
err_type <- c()
type <- c()
colors <- c()

name <- "conv-pangolin-sp-CFL0.7.dat"
dat_pango <- read.table(paste(folder_pango, name, sep=""))
append_data(dat_pango, dat_pango$V2, "corr", 1)

name <- "conv-pangolin-sp-nocorr-CFL0.7.dat"
dat_pango <- read.table(paste(folder_pango, name, sep=""))
append_data(dat_pango, dat_pango$V2, "nocorr", 3)


data <- data.frame(nb_cells, err, err_type, corr=type, colors)
data$corr <- factor(data$corr)
data$colors <- factor(data$colors) 

p <- ggplot(data, aes(x=nb_cells, y=err, colour=err_type)) 
p <- p + geom_line() + geom_point(aes(shape=corr), size=5)

# Legend
p <- p + scale_shape_manual(name="Wind correction", values=c(25,20), 
labels=c("with", "without"))
# Hide symbol in colour legend
p <- p + guides(colour = guide_legend(override.aes = list(size = 1)))
p <- p + scale_colour_discrete(name="Error", labels= 
    c(expression(l[2]), expression(l[oo]))) 

p <- p + xlab("nb cells") + ylab("error")
p <- p + theme(axis.title = element_text(size=16))
p <- p + ggtitle("Winds correction impact (sp, 0.28x0.18°)")
p <- p + coord_cartesian(ylim = c(1,2e-5)) 


# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(1e-1,1e-2,1e-3,1e-4,1e-5))
p <- p + scale_x_continuous(trans=log2_trans(),
                            breaks=c(1e3, 1e4, 1e5,1e6))
print(p)
dev.off()
