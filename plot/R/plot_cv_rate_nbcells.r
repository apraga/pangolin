# Compare cv rates vs nb cells for different models

library(ggplot2)
library(scales) # for log

source("nb_cells.r")

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
  err_type <<- c(err_type, rep("l[2]", size), rep("l[oo]", size))
  model <<- c(model, rep(nmodel, 2*size))
  type <<- c(type, rep(ntype, 2*size))
}

#-------------------------------------------------------------------------------

pdf("cv_rate_compare.pdf", width=10,height=5)
folder_pango <- "/wkdir/pae2/praga/parallel/cv_rate/"
folder <- "/wkdir/pae2/praga/lauritzen/conv/data/"
labels <- c("CLAW CFL0.95", "FARSIGHT CFL1.0", "SLFV-ML CFL0.8", 
    "UCISOM CFL0.8", "CAM-FV CFL0.2")

names <- c(
    "conv-CLAW-sp-gs-CFL0.95.dat",
    "conv-FARSIGHT-sp-gs-CFL1.0.dat",
    "conv-SLFV-ML-sp-gs-CFL0.8.dat",
    "conv-UCISOM-sp-gs-CFL0.8.dat",
    "conv-CAM-FV-sp-gs-CFL0.2.dat")

nb_cells <- c()
err <- c()
err_type <- c()
model <- c()
type <- c()

# Plot both unlimited and limited
for (ttype in c("sp", "un")) {
  name <- "conv-pangolin-sp-CFL0.96.dat"
  name <- gsub("sp", ttype, name)
  name <- paste(folder_pango, name, sep="")

  if (!file.exists(name)) {  next}
  dat_pango <- read.table(name)
  append_data(dat_pango, TRUE, dat_pango$V2, "PANGOLIN CFL0.95", ttype)

  for (i in seq(length(names))) {
    cur <- gsub("sp", ttype, names[i])
    cur <- paste(folder, cur, sep="")
    # Skip missing data
    if (!file.exists(cur)) { 
      print(cur)
        next 
    }
    dat <- read.table(cur)

    if (grepl("CLAW", names[i])) { n <- claw(dat$V1) }
    else if (grepl("FARSIGHT", names[i])) { n <- farsight(dat$V1) }
    else if (grepl("SLFV", names[i])) { 
    n <- c()
    for (l in seq(length(dat$V1))) {n[l] <- slfv(dat$V1[l]) }
    }
    else { n <- latlon(dat$V1)}

    append_data(dat, FALSE, n, labels[i], ttype)
  }
}
data <- data.frame(nb_cells, err, err_type, Model=model, type, parse=TRUE)

p <- ggplot(data, aes(x=nb_cells, y=err, colour=Model)) + geom_line(aes(linetype=type))
p <- p + geom_point(aes(shape=Model))
# We need facet_grid to parse the labels as mathematical expressions (does not
# work with facet_wrap)
p <- p + facet_grid(~err_type, labeller=label_parsed)#, ncol=2)  
# Custom facet labels
p <- p + theme(strip.text.x = element_text(size=12, face="bold"))

p <- p + xlab("Nb cells") + ylab("Error")
#p <- p + theme(axis.title = element_text(size=18))

# Custom Legend labels
p <- p + scale_linetype_discrete(name="Limiter", labels=c("shape-preserving",
"unlimited"))

# Log scale
p <- p + scale_y_continuous(trans=log10_trans(), labels=scientific_format(),
                            breaks=c(1e-1,1e-2,1e-3,1e-4,1e-5))
p <- p + scale_x_continuous(trans=log2_trans(),
                            breaks=c(1e3, 1e4, 1e5,1e6))
print(p)
dev.off()
