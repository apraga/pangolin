# Compare cv rates vs nb cells for hourdin test (2 config)

library(ggplot2)
library(scales) # for log

#-------------------------------------------------------------------------------
# Nb cells from resolution
pango <- function(lon) { 
  res <- floor(120/lon)
  res <- floor(0.5*(res+1))
  res <- 6*res^2
  return(res) 
}

# Append to global vectors
append_data <- function(dat, is_pango, nbcells, nmodel, ntype) {
  nb_cells <<- c(nb_cells, nbcells, nbcells)
  size <- length(dat$V2)
  # Pangolin has one more column
  err <<- c(err, dat$V3, dat$V4)
  err_type <<- c(err_type, rep("l2", size), rep("loo", size))
  model <<- c(model, rep(nmodel, 2*size))
  type <<- c(type, rep(ntype, 2*size))
}

#-------------------------------------------------------------------------------

pdf("cv_rate_hourdin.pdf", width=10,height=7)
folder_pango <- "/wkdir/pae2/praga/parallel/cv_rate/"

nb_cells <- c()
err <- c()
err_type <- c()
model <- c()
type <- c()

# Plot both unlimited and limited
for (ttype in c("sp", "hourdin_rot", "hourdin_norot")){#, "un")) {
  name <- "conv-pangolin-sp-CFL0.7.dat"
  name <- gsub("sp", ttype, name)

  dat_pango <- read.table(paste(folder_pango, name, sep=""))
  append_data(dat_pango, TRUE, dat_pango$V2, "PANGOLIN CFL0.7", ttype)

}
data <- data.frame(nb_cells, err, err_type, model, type)

levels(data$type) <- c("Hourdin", "Hourdin rotated", "Gaussian")
levels(data$model) <- c("Hourdin", "Hourdin rotated", "Gaussian")

p <- ggplot(data, aes(x=nb_cells, y=err, colour=type)) + geom_line(size=0.8)
p <- p + facet_wrap(~err_type, ncol=2)  
p <- p + xlab("nb cells") + ylab("error")
#p <- p + theme(axis.title = element_text(size=18))

# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(1e-1,1e-2,1e-3,1e-4,1e-5))
p <- p + scale_x_continuous(trans=log2_trans(),
                            breaks=c(1e3, 1e4, 1e5,1e6))
p <- p + ggtitle("Pole impact on convergence rate (Hourdin, one period)")
print(p)
dev.off()
