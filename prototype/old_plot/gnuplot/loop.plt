#set output 'tmp/'.N.'.gif' 
splot 'tmp/'.N.'.data' u 1:2:3:4 palette pt 5 ps 0.5
N = N +1
if (N < 2500) reread

