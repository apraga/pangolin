# Compare cv rates vs nb cells for different models
# This version plots for limited scheme, with hourdin comparison

library(ggplot2)
library(scales) # for log

# Append to global vectors
append_data <- function(dat, is_pango, nbcells, nmodel, ntype) {
  res <<- c(res, nbcells, nbcells)
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
labels <- c("UCISOM CFL1.0", "FARSIGHT CFL1.0", "CLAW", "CAM-FV",
            "SLFV-ML")


names <- c( "conv-UCISOM-sp-gs-CFL0.8.dat",
           "conv-FARSIGHT-sp-gs-CFL1.0.dat",
           "conv-CLAW-sp-gs-CFL0.95.dat",
           "conv-CAM-FV-sp-gs-CFL0.2.dat",
           "conv-SLFV-ML-sp-gs-CFL0.8.dat")


res <- c()
err <- c()
err_type <- c()
model <- c()
type <- c()

# Plot both unlimited and limited
for (ttype in c("sp", "un")) {
  name <- "conv-pangolin-sp-CFL0.7.dat"
  name <- gsub("sp", ttype, name)

  dat_pango <- read.table(paste(folder_pango, name, sep=""))
  append_data(dat_pango, TRUE, dat_pango$V1, "PANGOLIN CFL0.7", ttype)

  for (i in seq(length(names))) {
    #cur <- gsub("sp", ttype, names[i])
    cur <- names[i]
    cur <- paste(folder, cur, sep="")
    # Skip missing data
    if (!file.exists(cur)) { 
    print(cur)
    next }
    dat <- read.table(cur)

    append_data(dat, FALSE, dat$V1, labels[i], ttype)
  }
}
data <- data.frame(res, err, err_type, model, type)
p <- ggplot(data, aes(x=res, y=err, colour=model)) + geom_line(aes(linetype=type))
p <- p + facet_wrap(~err_type, ncol=2)  
p <- p + xlab("nb cells") + ylab("error")
#p <- p + theme(axis.title = element_text(size=18))
p <- p + ggtitle("Sp convergence rates")

# Custom reverse log
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
        log_breaks(base = base), 
        domain = c(1e-100, Inf))
}
# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(10, 1, 1e-1,1e-2,1e-3,1e-4,1e-5))
p <- p + scale_x_continuous(trans=reverselog_trans(2),
                            breaks=c(3, 1.5, 0.75, 0.375,0.1875))
p <- p + coord_cartesian(ylim = c(10e-5, 10))                            
#p <- p + scale_x_reverse()
print(p)
dev.off()
