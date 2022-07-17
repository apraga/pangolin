# Nb cells from resolution
latlon <- function(lon) { return(floor(360/lon)*floor(180/lon)) }

homme <- function(lon) { return(6*(90/(3*lon))**2) }

farsight <- function(lon) { return(6*(90/lon)**2) }

pango <- function(lon) { return(6*floor(0.5*(floor(120./lon)+1))**2)}

# Nb lat from nb cells
pango_nb_lat <- function(nb) {
  res <- floor(sqrt(nb/6.))
  # Find closest
  if (abs(6*(res+1)**2 - nb) < abs(6*res**2-nb)) { res <- res+1} 
  return(res)
}

lon <- c(3, 1.5, 0.75, 0.375, 0.1875)
mpas <- c(21506,86018,107522)
n_pango <- pango(lon)
n_homme <- homme(lon)
n_farsight <- farsight(lon)
#nb_lat <- c()
#for (i in lon) {
#  nb <- latlon(i)
#  nb_lat <- c(nb_lat, pango_nb_lat(nb))
#}
#print(sprintf("nb lat for latlon"))
#print(nb_lat)
#
#nb_lat <- c()
#for (i in lon) {
#  nb <- homme(i)
#  nb_lat <- c(nb_lat, pango_nb_lat(nb))
#}
#print(sprintf("nb lat for homme"))
#print(nb_lat)
#
#
#nb_lat <- c()
#for (i in lon) {
#  nb <- farsight(i)
#  nb_lat <- c(nb_lat, pango_nb_lat(nb))
#}
#print(sprintf("nb lat for farsight"))
#print(nb_lat)
#
#nb_lat <- c()
#nb_cells <- c(21506,86018,107522)
#for (nb in nb_cells) {
#  nb_lat <- c(nb_lat, pango_nb_lat(nb))
#}
#print(sprintf("nb lat for MPAS"))
#print(nb_lat)

#red <- 0.8*latlon
#patch <- 2*(90/lon)**2
#cubed2<- 6*(90/lon)**2
#voronoi<-c(0,21506,86018,0)
#gauss<-0.5*(360/lon)**2
#icos<-10*(180/(5*lon))**2+2
#voronoi2<-c(0,10000,20000,0)
