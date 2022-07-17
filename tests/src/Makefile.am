#-------------------------------------------------------------------------------
# Configuration for tests/sequential directory
#-------------------------------------------------------------------------------
if HAS_HDF5
bin_PROGRAMS = tests_sequential tests_parallel compare_hdf5
else
bin_PROGRAMS = tests_sequential tests_parallel 
endif

tests_sequential_SOURCES = tests_borders.f90 \
tests_cells.f90 \
tests_partitions.f90 \
tests_sequential.f90 \
tests_output.f90 

# Must also clean mod files
CLEANFILES = *.mod *.o

tests_borders.f90: tests_output.o
tests_cells.f90: tests_output.o
tests_partitions.f90: tests_output.o

#-------------------------------------------------------------------------------
# Configuration for tests/parallel directory
#-------------------------------------------------------------------------------
tests_parallel_SOURCES = tests_advection_2d.F90 \
tests_parallel.f90

# Compilation dependencies : include folder where the .mod files live
DEP_DIR=$(top_builddir)/src
FCFLAGS+=-I$(DEP_DIR)

if HAS_HDF5
compare_hdf5_SOURCES = compare_hdf5.c

compare_hdf5: compare_hdf5.o
	${H5CC} -o $@ $<

compare_hdf5.o: compare_hdf5.c
	${H5CC} -c $<
endif


# Linking dependencies : all objects file exect for the main program
src_simulation_OBJ=$(wildcard $(DEP_DIR)/*.o)
TEST_OBJ=$(filter-out $(DEP_DIR)/pangolin.o, $(src_simulation_OBJ))
tests_parallel_LDADD = $(TEST_OBJ)
tests_sequential_LDADD = $(TEST_OBJ)
compare_hdf5_LDADD = $(HDF5_LIBS)
