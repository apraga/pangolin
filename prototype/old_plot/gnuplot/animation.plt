N = 0
set pm3d
set view map
set terminal jpeg
#load 'loop.plt'
set output  for [i=1:2] 'pictures/pic'.i..jpg' 
splot 'tmp/'.i.'.data' u 1:2:3:4 palette pt 5 ps 0.5


