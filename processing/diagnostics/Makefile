SRC_DIAG=diag.f90 io.f90 pango_diagnostics.f90
SRC_EQ=extract_equator.f90
OBJ_DIAG=$(SRC_DIAG:.f90=.o)
OBJ_EQ=$(SRC_EQ:.f90=.o)

FC=gfortran
FCFLAGS=-g -fbacktrace -fcheck=all
H5FC=h5fc
H5LDFLAGS=-shlib

all: pango_diagnostics extract_equator

pango_diagnostics: ${OBJ_DIAG}
	${H5FC} ${H5LDFLAGS} -o $@  $^

extract_equator: ${OBJ_EQ}
	${FC} -o $@  $^

%.o: %.f90
	${FC} -c ${FCFLAGS} $<

io.o: io.f90
	${H5FC} -c $<

clean:
	rm -f ${OBJ_DIAG} ${OBJ_EQ}

exec:
	./pango_diagnostics

