# Create difference files for plotting convergence (see plot_mse.r)

dir <- "/wkdir/pae2/praga/convergence"

# Col 1 : nb lat, Col 2 : nb iterations
# Change here (and use zeros padding)
names <- c('025', '030', '045', '090', '135', '180', '225', '270', '360', '450',
           '675',
           '220', '263', '392', '777', '1163', '1549', '1935', '2320', '3092', 
           '3864', '5792')
names <- matrix(names, ncol=2)

# Subset
nstart <- 1
nend <- dim(names)[1]

# Find file name for a given iteration
# Assume a certain format
findname <- function(i, names, ref) {
  if (ref == 1) {
    cur <- paste(dir, "/ratio_anal_", names[i,1], "_", names[i,2], ".dat", sep="")
  }
  else {
    cur <- paste(dir, "/sequential/ratio_seq_", names[i,1], "_", names[i,2], ".dat", sep="")
  }
  return (cur)
}

# Create lists of vector, append with diff <- c(diff, list(myvar))
# Here we only write to an output file
for (i in seq(nstart, nend)) {
  ref <- read.table(findname(i, names, 1))
  q <- read.table(findname(i, names, 0))
  diff <- ref[1] - q[1]
  fname <- paste(dir, "/diff_seq_", names[i,1], "_", names[i,2], ".dat", sep="")
  write(unlist(diff), file=fname,sep="\n") 
}
