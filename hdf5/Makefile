#FC=/usr/lib64/openmpi/bin/h5pfc 
FC=h5pfc 

pango_write: pango_write.o
	${FC} -g -shlib -o $@ $<

pango_write.o: pango_write.f90
	${FC} -g -c $<
