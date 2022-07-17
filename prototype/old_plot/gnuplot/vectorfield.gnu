scale = 500
plot '../data/data_vectorfield.dat' u 1:2:($3*scale):(-$4*scale) w vectors
