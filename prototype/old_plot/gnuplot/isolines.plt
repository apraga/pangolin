#set pm3d
#set view map
set xrange [0:2*pi]
set yrange [0:pi/2]
splot cos(0.5*x)**2*cos(y)**2 with lines
#splot 'tmp/'.1.'.data' u 1:2:3:4 palette pt 5 ps 0.5



#set pm3d
#set parametric
#set urange [0:2*pi]
#set vrange [0:pi/2]
## Parametric functions for the sphere
#fx(v,u) = cos(v)*cos(u)
#fy(v,u) = cos(v)*sin(u)
#fz(v)   = sin(v)
#color(u,v) = cos(0.5*u)*cos(v)*2 
#
##splot fx(v,u),fy(v,u),color(u,v)
##splot fx(v,u),fy(v,u),fz(v)
#splot fx(v,u),fy(v,u),color(u,v)
