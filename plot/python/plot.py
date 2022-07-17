import sys 
import lib_plot # personal library
import getopt


###############################################################################
# Wrapper for plotting function
###############################################################################
try:
  argv = sys.argv[1:]
  opts, args = getopt.getopt(argv,"p:if:ew")
except getopt.GetoptError:
  print('plot.py -p <projection> [-i] -f <inputfile> [-e] [-w]')
  print('-i: activates interpolation')
  print('-p: miller, ortho, npolor, spolar, zoom_ortho')
  print('-e: for an error plot')
  print('-w: for a winds plot')
  sys.exit(2)

proj = None
interpol = None
fname = None
iserr = None
iswinds = None
for opt, arg in opts:
  if opt in ("-h"):
    print('plot.py -p <projection> -i <interpolation> -f <inputfile> -e')
    sys.exit(2)
  elif opt in ("-p"):
    proj = arg
  elif opt in ("-i"):
    interpol = True
  elif opt in ("-f"):
    fname = arg
  elif opt in ("-e"):
    print("Plotting error")
    iserr = arg
  elif opt in ("-w"):
    print("Plotting winds (needs zonal as input)")
    iswinds = True

#title = "Numerical - analytical (rotated snail test, T, 0.56x0.38$^\circ$)" 
title = "Numerical - analytical (gaussian hills, T, 0.56x0.38$^\circ$)" 
#title = "Numerical (rotated snail test, T, 0.56x0.38$^\circ$)" 
name_fig = 'switched.png'#pango-seq-err-gs-160lat-T.png'

if fname is None:
  print("Must specify input file")
  sys.exit(2)

lib_plot.plot_data(fname, proj, interpol, iserr, iswinds, title, name_fig)
