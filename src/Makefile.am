#-------------------------------------------------------------------------------
# Configuration for source files, recursive.
#-------------------------------------------------------------------------------
#
bin_PROGRAMS = pangolin
noinst_PROGRAMS = pangolin

#FCFLAGS+=-I/space/praga/softwares/mpe/include 
#pangolin_LDADD=-L/space/praga/softwares/mpe/lib  -lmpe_f2cmpi -llmpe -lmpe -lpmpi -lmpi

pangolin_LDADD=@PROFILING_LIBS@ 

# Must also clean mod files
CLEANFILES = *.mod 

dist_pangolin_SOURCES = advection_2d.F90 \
analytical_partitioning.F90 \
analytical_data.inc \
band_grid_class.F90 \
configuration_class.F90 \
data_check.F90 \
date.c \
debug.f90 \
global_grid_class.f90 \
io.inc \
io_ascii.inc \
io_hdf5.inc \
list_tracers_class.F90 \
message.F90 \
parameters.F90 \
partition_class.F90 \
partitioning_class.F90 \
profiling.F90 \
pangolin.F90 \
pangolin_run.F90

# Additional dependencies
pangolin.F90 : pangolin_run.o
partitioning_class.F90: partition_class.o global_grid_class.o \
  band_grid_class.o analytical_partitioning.o profiling.o
partition_class.F90: band_grid_class.o global_grid_class.o \
	analytical_partitioning.o analytical_data.inc 
band_grid_class.o: global_grid_class.o analytical_partitioning.o \
  message.o parameters.o list_tracers_class.o
configuration_class.F90: data_check.o
advection_2d.F90: partition_class.o partitioning_class.o \
  parameters.o
global_grid_class.f90: configuration_class.o
data_check.F90: message.o parameters.o
pangolin.F90 : pangolin_run.o
partitioning_class.F90: partition_class.o global_grid_class.o \
  band_grid_class.o analytical_partitioning.o
partition_class.F90: band_grid_class.o global_grid_class.o \
	analytical_partitioning.o analytical_data.inc 
band_grid_class.o: global_grid_class.o analytical_partitioning.o \
  message.o parameters.o
configuration_class.F90: data_check.o
advection_2d.F90: partition_class.o partitioning_class.o \
  parameters.o
global_grid_class.f90: configuration_class.o
data_check.F90: message.o
advection_2d.F90: profiling.o

advection_2d.o: date.o
configuration_class.o:: date.o

# Generate ctags
ctags:
	ctags --langmap=Fortran:.inc.f90.F90 *.{f90,F90,inc}
