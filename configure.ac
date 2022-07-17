AC_INIT([pangolin], [0.1], [praga@cerfacs.fr])

# Auxiliary directory
AC_CONFIG_AUX_DIR([build-aux])

# Foreign is for non-GNU syntax
# Also, parallel tests
m4_version_prereq([2.66],[
  AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests])
  ],[
  AM_INIT_AUTOMAKE([subdir-objects foreign])
])


list_compil='mpif90 mpiifort ifort'

dnl If the user specifies it, search in path a intel, pgf or gfortran MPI compiler
AC_ARG_WITH([compiler],
            [AS_HELP_STRING([--with-compiler=COMPIL], [set compiler with MPI
             support (disabled by default)])],
            [search_mpi=$withval])

dnl Enable debug flags
 AC_ARG_ENABLE([debug],
             AS_HELP_STRING([--enable-debug], [enable debug mode (traceback etc.)])
             [], has_debug=1)

dnl Enable openmp
 AC_ARG_ENABLE([openmp],
             AS_HELP_STRING([--enable-openmp], [activates OpenMP])
             [], has_openmp=1)

dnl Enable profiling
 AC_ARG_ENABLE([profiling],
             AS_HELP_STRING([--enable-profiling], [activates profiling. Intel by 
                             default. If MPE_DIR is set, will use MPE instead.])
             [], has_profiling=1)

 dnl Enable optim
 AC_ARG_ENABLE([optim],
              AS_HELP_STRING([--enable-optim], [activates
                              optimization.])
              [], has_optim=1)
 dnl Enable ascii and disable netcdf
 AC_ARG_ENABLE([ascii],
              AS_HELP_STRING([--enable-ascii], [activates
                              sequential ASCII I/O and disable netCDF I/O.])
              [], has_ascii=1)



dnl MPI_WRAPPER overrides anything else
if test x"$MPI_DIR" = x; then
  if test x"$search_mpi" = x; then
    AC_MSG_NOTICE([No compiler specified, default is gfortran])
    search_mpi="gfortran"
fi

  AC_MSG_CHECKING(for $search_mpi compiler with MPI support)
  GET_MPI_COMPILER($list_compil, $search_mpi)
else
  MPI_WRAPPER=mpif90
  MPI_RUN=mpirun
  AC_MSG_CHECKING(for $search_mpi compiler with MPI support)
fi
 
dnl Set the FC variable to the previously found value. If the search failed,
dnl search for a default executable
if test x"$MPI_WRAPPER" != x; then
  if test x"$MPI_DIR" != x; then
    FC=$MPI_DIR/$MPI_WRAPPER
    CC=$MPI_DIR/mpicc
  else
    FC=$MPI_WRAPPER
    CC=mpicc
  fi
  AC_MSG_RESULT($MPI_WRAPPER ($MPI_DIR))
else
  AC_MSG_RESULT(failed)
fi
 
dnl No compiler specified or found, try a default script
if test x"$MPI_WRAPPER" = x; then
  AC_MSG_NOTICE([No compiler specified or found in the PATH, will use default
                 MPI wrapper])
  dnl If not, we try our last chance : the script listed in $list_compil
  dnl directly
  for run_test in $list_compil; do
    AC_MSG_CHECKING([for default MPI script $run_test])
    if_exists=`echo "$run_test" | grep "not found"`
    if test x"$if_exists" = x; then
      AC_MSG_RESULT([exists])
      FC=$run_test
      last_chance=1
      break
    else
      AC_MSG_RESULT([does not exist])
    fi
  done
fi

if test x"$FC" = x; then
  AC_MSG_ERROR([could not find Fortran compiler for MPI])
fi

dnl Cannot use AC_LANG([Fortran]) 
AC_PROG_FC 
AC_PROG_CC 
AM_PROG_CC_C_O

dnl Temporary variable if MPI_WRAPPER is overriden
compiler=$search_mpi
if test x$compiler == x; then
  GET_MPI_COMPILER_VERSION($FC, compiler)
fi

dnl Flags for compiler. Adds eventual debug and openmp flags.
dnl As Fortran support is not great, must do it by hand
case $compiler in
  "pgfortran") 
    PREPROC_FLAGS="-DWITH_PGF90"

    if test x"$has_debug" != x; then
      DEBUG_FLAGS="-C -g -O0 -Mprof=mpich2,lines" #for profiling
    fi
    if test x"$has_openmp" != x; then
      OPENMP_FLAGS="-mp"
    fi
  ;;
  "ifort")
    PREPROC_FLAGS="-DWITH_IFORT"
    if test x"$has_debug" != x; then
      DEBUG_FLAGS="-traceback -check all -fp-stack-check"
    fi
    if test x"$has_openmp" != x; then
      OPENMP_FLAGS="-openmp"
    fi
  ;;
  "gfortran")
    if test x"$has_debug" != x; then
      DEBUG_FLAGS="-Wall -Wextra -fimplicit-none -fbounds-check -fbacktrace "
    fi
    if test x"$has_openmp" != x; then
      OPENMP_FLAGS="-fopenmp"
    fi
    ;;
  *) ;;
esac

dnl Add custom debug flags and -g
AC_MSG_CHECKING([for debug mode])
if test x"$has_debug" != x; then
  AC_MSG_RESULT(yes)
  DEBUG_FLAGS+=" -DDEBUG -g"
else
  AC_MSG_RESULT(no)
fi

dnl Cannot use AC_OPENMP
if test x"$has_openmp" != x; then
AC_MSG_CHECKING([for OpenMP flags])
AC_COMPILE_IFELSE([
            program main
            end],[echo "$OPENMP_FLAGS"],[echo "failed"])
fi


dnl profiling
AC_MSG_CHECKING(for profiling mode)
if test x"$MPE_DIR" != x ; then
  AC_MSG_RESULT([yes, MPE])
  PROFIL_FLAGS=-I$MPE_DIR/include
  dnl Must be after the objects
  PROFILING_LIBS="-L$MPE_DIR/lib -lmpe_f2cmpi -llmpe -lmpe -lmpi -lpmpi"
  AC_SUBST(PROFILING_LIBS)
else 
  if [ test x"$has_profiling" != x && test x"$compiler" == x"ifort" ]; then
    PROFIL_FLAGS=-tcollect
    AC_MSG_RESULT([yes, Intel])
  else
    AC_MSG_RESULT([no, need Intel MPI or MPE])
  fi
fi

dnl Optimization
AC_MSG_CHECKING(for optimization (override profiling))
if [ test x"$has_optim" != x ]; then
  PROFIL_FLAGS=-O3
  AC_MSG_RESULT(yes)
 else
  AC_MSG_RESULT(no)
fi

AC_MSG_CHECKING(for I/O type)
if [ test x"$has_ascii" == x ]; then
   AC_MSG_RESULT(hdf5)
 else
  AC_MSG_RESULT(sequential ASCII)
fi

dnl dnl netCDF configuration
dnl if [ test x"$has_ascii" == x ]; then
dnl 
dnl   dnl Search for netCDF executables which will give us directories
dnl   AC_PATH_PROG([NF_CONFIG], [nf-config], [], [$NETCDF_FDIR])
dnl   AC_PATH_PROG([NC_CONFIG], [nc-config], [], [$NETCDF_CDIR])
dnl   if [ test x"$NC_CONFIG" == x ]; then
dnl     AC_MSG_FAILURE([did not found nc-config])
dnl   fi
dnl 
dnl   dnl Define flags for Fortran netcdf
dnl   dnl Old versions do not have fortran and c packaged separetely
dnl   NF_CONFIG2=$NF_CONFIG
dnl   if [ test x"$NF_CONFIG" == x ]; then
dnl     AC_MSG_NOTICE([C and Fortran version are not separated])
dnl     NF_CONFIG2=$NC_CONFIG
dnl   fi
dnl   AC_SUBST([NETCDF_FINC],[$( $NF_CONFIG --includedir)]) 
dnl   NETCDF_FCFLAGS=" -I$NETCDF_FINC -DNETCDF"
dnl   AC_SUBST([NETCDF_FLIBS],[$( $NF_CONFIG --flibs)]) 
dnl 
dnl   dnl Define flags for C netcdf
dnl   AC_SUBST([NETCDF_CINC],[$( $NC_CONFIG --includedir)]) 
dnl   NETCDF_CFLAGS=" -I$NETCDF_CINC"
dnl   AC_SUBST([NETCDF_CLIBS],[$( $NC_CONFIG --libs)]) 
dnl fi

AM_CONDITIONAL(HAS_HDF5, [test x"$has_ascii" == x])
if [ test x"$has_ascii" == x ]; then
  SET_HDF5()
fi


FCFLAGS="$PREPROC_FLAGS $DEBUG_FLAGS $OPENMP_FLAGS $PROFIL_FLAGS $HDF5_FLAGS"
LIBS="$HDF5_LIBS"
CFLAGS+=" $HDF5_FLAGS"

dnl Find executable for simulation and store it into MPI_RUN
AC_ARG_VAR(MPI_RUN,[starts an MPI job])
AC_ARG_VAR(MPI_DIR,[MPI directory])
AC_PATH_PROG(MPI_RUN, qsub)

dnl Define a conditional. It will serve in Makefile.am to create the executable
dnl script. We need to know if we can use qsub
AM_CONDITIONAL(IS_CLUSTER, [test x"$MPI_RUN" != x])

AC_MSG_CHECKING(for an executable for starting MPI jobs)
if test x"$MPI_RUN" = x; then
  dnl If default mpi compiler, there is no directory
  if test x"$MPI_DIR" = x; then
    MPI_RUN=mpirun
  else
    MPI_RUN=$MPI_DIR/mpirun
  fi
fi

if test x"$MPI_RUN" = x; then
  AC_MSG_ERROR([could not find script for MPI jobs])
fi
AC_MSG_RESULT(mpi_run ($MPI_DIR))

dnl if Test::More is too old, install cpanm and the local::lib module in ~/perl5
dnl Then we can use the proper version
dnl AC_PERL_MODULE_VERSION(Test::More 1.0, , wget -O- http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib)

dnl Default installation directory is here
AC_PREFIX_DEFAULT([`pwd`])

dnl Create needed directories is here
AC_MSG_CHECKING(for log and output directories)
AC_CONFIG_COMMANDS([log], [mkdir -p log])
AC_CONFIG_COMMANDS([output], [mkdir -p output])
AC_MSG_RESULT(ok)

AC_CONFIG_HEADERS([config.h])
dnl Recursive for src and tests
AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
                 tests/Makefile
                 tests/src/Makefile
                 generate/Makefile
                 ])

AC_OUTPUT
