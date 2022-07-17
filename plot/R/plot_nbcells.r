# Plot total nb of interior cells from a file containing different cell numbers.
# File format : interior/halo cells zonal flux, interior/halo gradient for
# meridional, interior/hal for meridional flux, interior/halo for meridional
# slope, total interior cells

dir <- "/wkdir/pae2/praga/nb_cells/90lat"

nb_procs <- seq(3,63,3)

nb_cells <- rep(1, length(length(nb_procs)))

i <- 1
for (n in nb_procs) {
    cur <- paste(dir, "/nb_cells", n, sep="")
    nb <- read.table(cur)
    nb_cells[i] <- max(nb[,10])
    i <- i+1
}

postscript("nbcells.ps",width=10,height=10,pointsize=15)
plot(nb_procs, nb_cells, type="o",xlab="cores", ylab="interior cells")

# Theoritical
par(new=TRUE)
x <- seq(3,63)
y <- 90^2 / (x/6)
extrem=c(min(nb_cells), max(nb_cells))
plot(x, y,ylim=extrem, type="l", lty=2, ylab="", xlab="")
axis(1, at=c(3),labels="3")
title("Domain size (1x0.67° at the Equator)")
legend("topright", c("Numerical", "Theorical"), lty=c(1,2))
dev.off()
