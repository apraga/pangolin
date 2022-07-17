# Compute the maximal timestep according to the CFL condition
# (for zonal avection)
#set terminal postscript eps
#set output  'condition_cfl_zonal.eps' 

set xrange [1:180] noreverse
set yrange [1:1000] 
set samples 180
coef = pi/180
f(x) = x > 90 ? 2*(181-x)-1 : 2*x-1
plot 180*(2*pi)/(cos((90.5-x)*coef)*3*f(x))
#plot cos((90.5-x)*coef)*3*f(x)/(2*pi)

