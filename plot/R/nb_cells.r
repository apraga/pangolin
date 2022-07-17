# Compute the number of cells for some models given a resolution

latlon <- function(lon) { return(floor(360./lon)*floor(180./lon)) }

#homme <- function(lon) { return(6*(90/(3*lon))**2) }

farsight <- function(lon) { return(6*(90./lon)**2) }

claw <- function(lon) { return(2*(90./lon)**2) }

slfv <- function(lon) {
 val <- c(5762,23042,92162,368642)
 if (lon == 3.) { i = 1}
 else if (lon == 1.5) { i = 2}
 else if (lon == 0.75) { i = 3}
 else if (lon == 0.375) { i = 4}
 else { return(0)}
 return(val[i])
}

pango <- function(lon) { 
  res <- ceiling(120/lon)
  # Upper integer for more precision
  res <- ceiling(0.5*(res+1))
  res <- 6*res^2
  return(res) 
}

# Find nb lat2 from nb cells
pango_nb <- function(nb) { 
  res <- ceiling(sqrt(nb/6))
  return(res) 
}
#mpas <- function(lon) { 
#  res <- pango(lon)
#  return(2.24*res) 
#}

# Snippet for finding nb_lat2 for Pangolin to be as close as possible of the
# number of cells of each model
#res <- c(3.0000 1.5000 0.7500 0.3750 0.1875) 
res <- 0.75
n <- latlon(res)
#print(sprintf("for latlon %f", pango_nb(n)))
n <- claw(res)
#print(sprintf("for claw %f", pango_nb(n)))
n <- slfv(res)
#print(sprintf("for slfv %f", pango_nb(n)))
n <- farsight(res)
#print(sprintf("for farsight %f", pango_nb(n)))
