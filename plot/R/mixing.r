# Contains function for plotting mixing

append_data <- function(ttype, name_in, mmodel, ffolder, sskip) {
  name <- gsub("sp", ttype, name_in)
  name <- paste(ffolder, name, sep="")
  if (!file.exists(name)) { return() }
  q <- read.table(name, skip=sskip)
  khi <<- c(khi, q$V1)
  xi <<- c(xi, q$V2)
  model <<- c(model, rep(mmodel, length(q$V1)))
  if (ttype == "sp") { tname <- "shape preserving" }
  else { tname <- "unlimited" }
  type <<- c(type, rep(tname, length(q$V1)))

}

create_data <- function() {
  # Plot both unlimited and limited
  for (ttype in c("sp", "un")) {
    for (i in seq(length(names))) {
      sskip <- 0
      # Special header for farsight
      if (grepl("farsight", names[i])) { sskip <- 10}
      print(sprintf("skid %d names %s\n", sskip, names[i]))
      append_data(ttype, names[i], models[i], folder, 10)
    }
  }
}

create_data_pango <- function() {
  n <- length(names)
  for (ttype in c("sp", "un")) {
    append_data(ttype, names_pango, "Pangolin CFL0.95", folder_pango, 0)
  }
}
