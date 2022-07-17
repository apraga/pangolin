# Write grid data on VTK format
import numpy

f = open('reduced_grid.vtk','w')
# Header
f.write("# vtk DataFile Version 1.0 \n")
f.write("Reduced grid\n")
f.write("ASCII\n")

nb_lat = 2
nb_points = 3*nb_lat*nb_lat
f.write("DATASET UNSTRUCTURED_GRID\n")

f.write("POINTS "+str(nb_points)+" float\n")
f.close()

lat = 0.
lon = 0.
x = numpy.cos(lon)*numpy.sin(lat)
y = numpy.sin(lon)*numpy.sin(lat)
z = numpy.cos(lat)

s = "%f %f %f \n" % (x,y,z)
f.write(s)

i = 1
nb_lon = 3*(2*i-1)
dlat = 360./nb_lon
#    s = "%f %f %d %d \n" % (middle,1/step,1,i)
#    #f.write(s)
#split_single(3.0,list_middle,f)
#
#for k in range(2,91):
# split(k,f)


