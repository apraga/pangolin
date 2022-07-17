DIR=/space/praga/softwares/scotch_5.1.12
FCFLAGS=-I${DIR}/include
LDFLAGS=-L${DIR}/lib -lscotch -lscotcherr -lm

all: partitioning

partitioning: grid.o output.o partitioning.o 
	gfortran -o $@ $^ ${LDFLAGS}

%.o: %.f90
	gfortran ${FCFLAGS} -c $<

clean: 
	rm -f *.o *.mod
