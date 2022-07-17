# Compile these folders
SUBDIRS = src tests generate

# Contains customs macros
ACLOCAL_AMFLAGS = -I m4

# There we add misc files needed for the correct execution : batch and
# configuration templates
EXTRA_DIST = pbs_run.batch tests/pbs_run.batch config.template \
	     simulation_run.pl SimulationRun.pm config \
	     tests/functional.t tests/parallel.t tests/sequential.t \
	     tests/TestsRun.pm tests/CompareFiles.pm tests/tests_run.pl


#-------------------------------------------------------------------------------
# This is done as a custom hook for bin_PROGRAMS would put it into a bin/ 
# directory.
#-------------------------------------------------------------------------------
install-exec-hook: update_scripts

# Create temporaray folder with several configuration files (one for each
# situation). We consider the max
TMP=tmp
NB_LAT2=90

# Non-blocking tests
QSUB=qsub

#-------------------------------------------------------------------------------
# Update the scripts with the correct mpi executable
#-------------------------------------------------------------------------------
BATCH=pbs_run.batch
EXE_SCRIPT=SimulationRun.pm

update_scripts:
if IS_CLUSTER
LOCRUN=$(QSUB)
LOCBATCH=$(BATCH)
else
LOCRUN=$(MPI_RUN)
endif
	for cur in $(EXE_SCRIPT); do \
	  echo $$cur;\
	  sed -i 's:$$mpirun=.*:$$mpirun="$(LOCRUN)";:g' $$cur ;\
	  sed -i 's:$$batch=.*:$$batch="$(LOCBATCH)";:g' $$cur ;\
	  sed -i 's:$$mpidir=.*:$$mpidir="$(MPI_DIR)";:g' $$cur ;\
	chmod +x $$cur ;\
	done

# Force to rebuild
.PHONY: update_scripts

#-------------------------------------------------------------------------------
# Below are our custom functions
#-------------------------------------------------------------------------------

# Cleaning
cleanlog:
	rm -f log/* pangolin.o* pangolin.log*

cleantrace:
	rm -f simulation.*

# Send code to neptune
sendneptune:
	@rsync --include-from=rsync.filter -avz . praga@neptune:phd/pangolin-0.1/

# Send code to corail
sendcorail:
	@rsync --include-from=rsync.filter -avz . praga@corail:phd/pangolin-0.1/

# Receive simulation results from neptune
receiveneptune:
	@rsync --include="*time*" --exclude="*" -avz \
	  praga@neptune:phd/code/output/ output/ 

# Receive simulation results from corail
receivecorail:
	@rsync --include="*time*" --exclude="*" -avz \
	  praga@corail:phd/code/output/ output/ 

# Check no run has failed
checkrun:
	@if grep -l -i -e termin -e fail -e segmentation log/*; \
	  then echo "Some runs failed"; \
	else echo "Runs ok"; \
	fi
# Custom
valgrind:
	LD_PRELOAD=/usr/lib64/valgrind/libmpiwrap-amd64-linux.so /usr/lib64/openmpi/bin/mpirun -np 6 valgrind --log-file=valgrind.log -v src/simulation
