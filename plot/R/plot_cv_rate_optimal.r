library(ggplot2)

# Compute optimal cv rate by linear regression

# Specify the column containing the error
append_data <- function(name, col) {
  if (!file.exists(name)) { 
    print(name)
    return()
  }

  q <- read.table(name)

  for (k in c(col, col+1)) {
    # As Lauritzen, use a non-linear fit, with starting point
    # For Pangolin, use asymptotic regime data, as CLAW does (that is, the last
    # 3 data points)
    y <- log(q[[k]][4:6])
    x <- log(q[[1]][4:6])
    fit <- nls(y~p1*x+p2, start=list(p1=1.1, p2=1.1))

    # Extract slope
    coef <<- c(coef, summary(fit)$parameters[1,1])
    type <<- c(type, ttypes[j])
    # Custom spacing
    pos <<- c(pos, 3*i+j)
    # For several plots
    if (k == col) { err_type <<- c(err_type, "l[2]") }
    else { err_type <<- c(err_type, "l[oo]") }
  }
}

# Read data from txt file. We read line by line as the number of column differ
read_histogram_data <- function(input) {
  con  <- file(input, open = "r")

  i <- 1
  # Colums are separated by tab
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    vec <- (strsplit(oneLine, "\t"))

    # Skip vector of empty strings
    if (sum(nchar(vec[[1]])) == 0) { next}

    # Skip if model not in the list
    tmp <- grepl(vec[[1]][2], list_names)
    if (!any(tmp == TRUE)) { next}

    labels <<- c(labels, vec[[1]][2])
    print(vec[[1]][2])
    # Extract convergence rate to numerical et extract vector
    vec <- list(as.numeric(vec[[1]][5:8]))[[1]]
    # Remove negative values
    vec <- replace(vec, vec < 0, 0)
    
    # Data is ordered : k2 un, koo un, k2 sp , koo sp
    coef <<- c( coef, vec)
    type <<- c(type, "un", "un", "sp", "sp")
    err_type <<- c(err_type, "l[2]", "l[oo]", "l[2]", "l[oo]")
    pos <<- c(pos, 3*i+1, 3*i+1, 3*i+2, 3*i+2 )
    i <- i + 1
  } 
  close(con)
}
#-------------------------------------------------------------------------------

# Numerical convergence order
folder <- "/wkdir/pae2/praga/lauritzen/conv/data/"
folder_pango <- "/wkdir/pae2/praga/parallel/cv_rate/"
names <- c(
    "conv-CLAW-sp-gs-CFL0.95.dat",
    "conv-FARSIGHT-sp-gs-CFL1.0.dat",
    "conv-SLFV-ML-sp-gs-CFL0.8.dat",
    "conv-UCISOM-sp-gs-CFL0.8.dat",
    "conv-CAM-FV-sp-gs-CFL0.2.dat")
names_pango <- "conv-pangolin-sp-CFL1.0.dat"
labels <- c()
list_names <- c("CLAW-CN0.95", "FARSIGHT-CN1.0", "SLFV-ML-CN0.8", 
    "UCISOM-CN0.8", "CAM-FV-CN0.2")
 
#labels <- c("CLAW CFL0.95", "FARSIGHT CFL1.0", "SLFV-ML CFL0.8", 
#    "UCISOM CFL0.8", "CAM-FV CFL0.2", "PANGOLIN CFL1.0")
ttypes <- c("sp", "un")

coef <- c()
type <- c()
pos <- c()
err_type <- c()

# Histogram data is read from text file
input <- "/wkdir/pae2/praga/lauritzen/histogram/all.txt"
read_histogram_data(input)


i <- length(list_names)+1
# Pangolin data
for (j in seq(length(ttypes))) {
  name <- gsub("sp", ttypes[j], names_pango[1])
  name <- paste(folder_pango, name, sep="")
  # Offset here
  append_data(name, 3)
  labels <- c(labels, "PANGOLIN-CN0.95")
}

pdf("opt_cv_rate_compare.pdf", width=10,height=5)
data <- data.frame(pos, coef, type, err_type)
p <- ggplot(data, aes(x=pos, y=coef, fill=type))
# Set histogram side by side
p <- p + geom_histogram(position="dodge", stat="identity")

# Custom breaks
breaks <- 3*seq(length(labels))+1.5
p <- p + scale_x_continuous(labels=labels, breaks=breaks)
p <- p + facet_grid(~err_type, labeller=label_parsed)
# Custom facet labels
p <- p + theme(strip.text.x = element_text(size=12, face="bold"))


# rotate text
p <- p + theme(axis.text.x = element_text(angle = 30), 
    axis.title.x = element_blank())

# Customise legend
p <- p + scale_fill_discrete(name="Limiter", labels=c("shape-preserving",
"unlimited"))
p <- p + ylab("Optimal convergence rate")

print(p)
dev.off()
