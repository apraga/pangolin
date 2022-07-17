# SET_HDF5
# Search if HDF5 is availabe and define the proper variables and flags
#
AC_DEFUN([SET_HDF5],
[
  dnl hdf5 configuration
  dnl Search for HDF5 executables which will give us directories
  dnl search first HDF5_BINDIR, then the PATH
  if [ test x"$HDF5_BINDIR" == x ]; then
    HDF5_BINDIR=$HDF5_DIR/bin
  fi
  
  AC_PATH_PROG([H5FC], [h5pfc], [], [$HDF5_BINDIR$PATH_SEPARATOR$PATH])
  if [ test x"$H5FC" == x ]; then
    AC_PATH_PROG([H5FC], [h5fc], [], [$HDF5_BINDIR$PATH_SEPARATOR$PATH])
  fi
  
  if [ test x"$H5FC" == x ]; then
    AC_MSG_FAILURE([did not found h5fc or h5pfc. Is HDF5 really installed ?])
  fi
  AC_PATH_PROG([H5CC], [h5cc], [], [$HDF5_BINDIR$PATH_SEPARATOR$PATH])
  if [ test x"$H5CC" == x ]; then
    AC_PATH_PROG([H5CC], [h5pcc], [], [$HDF5_BINDIR$PATH_SEPARATOR$PATH])
  fi
  dnl Compile with shared libraries
  H5CC="$H5CC -shlib"

  if [ test x"$H5CC" == x ]; then
    AC_MSG_FAILURE([did not found h5cc or h5pcc. Is HDF5 really installed ?])
  fi
   
  dnl HDF5_DIR override the search
  if [ test x"$HDF5_DIR" != x ]; then
    HDF5_FLAGS+=" -I$HDF5_DIR/include"
    HDF5_CLAGS+=" -I$HDF5_DIR/include"
    HDF5_LIBS+=" -Wl,-rpath,$HDF5_DIR/lib -L$HDF5_DIR/lib"
    HDF5_LIBS+=" -lhdf5_fortran -lhdf5 -lhdf5_hl -lz"
  
  else
    dnl Extract compiler flags from h5fc or h5pfc. 
    dnl The rule is : "everything from -I to -L"
    AC_SUBST([TMP_FFLAGS],[$($H5FC -show -shlib | perl -lne 'print for m/(-I\/.*)-L/g')]) 
    AC_SUBST([TMP_CFLAGS],[$($H5CC -show -shlib | perl -lne 'print for m/(-I\/.*)-L/g')]) 
    dnl Extract linker flags from h5fc or h5pfc. 
    dnl The rule is : "everything after -L to the end"
    AC_SUBST([HDF5_LIBS],[$($H5FC -show -shlib | perl -lne 'print for m/(-L.*)/g')]) 
 
    dnl Allow the user to add its flags 
    HDF5_FLAGS+="$TMP_FFLAGS"
    HDF5_CFLAGS+="$TMP_FFLAGS"
  
  fi
  HDF5_FLAGS+=" -DHDF5"

  dnl Set compilation flags temporarily 
  FCFLAGS="$HDF5_FLAGS"
  CFLAGS="$HDF5_CFLAGS"
  LIBS="$HDF5_LIBS"

  dnl Fortran language with fortrang 90 extensions for the following tests
  AC_LANG_PUSH([Fortran])
  AC_FC_SRCEXT([f90])

  dnl Try to compile a program to check HDF5 support
  AC_MSG_CHECKING(whether we can compile and link with HDF5)
  AC_LINK_IFELSE(
      [AC_LANG_SOURCE([
      [program main
       use hdf5
       contains
       end program]])],
      [AC_MSG_RESULT([ok])], 
      [AC_MSG_FAILURE([Could not compile and link with HDF5])])
  
  dnl Check if HDF5 is compiled with parallel
  dnl Extract the list of included dirs, the word after each -I
  list_dirs=`echo $HDF5_FLAGS | perl -lne 'print for m/-I(\S+)(\s|$)/g'`
  dnl For each directory, search if HAS_PARALLEL 1 is in one header
  for cur in $list_dirs; do
    if test -f $cur/H5config.h; then
      tmp=`grep "HAVE_PARALLEL 1" $cur/H5config.h`
    elif test -f $cur/H5pubconf.h; then
      tmp=`grep "HAVE_PARALLEL 1" $cur/H5pubconf.h`
    fi
    PARALLEL_HDF5=$PARALLEL_HDF5$tmp
  done
  AC_MSG_CHECKING(if HDF5 is parallel)
  if [ test x"$PARALLEL_HDF5" != x]; then
    HDF5_FLAGS+=" -DPARALLEL_HDF5"
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_WARN([HDF5 is not parallel])
    AC_MSG_WARN([This version does not work with non-parallel HDF5 !])
  fi

  dnl Check if the Fortran module truly has parallel support
  if [ test x"$PARALLEL_HDF5" != x]; then
    AC_MSG_CHECKING(whether we can compile and link with parallel HDF5)
    AC_COMPILE_IFELSE(
      [AC_LANG_SOURCE(
        [program test
         use hdf5
         implicit none
         include 'mpif.h'
         character(7), parameter :: filename = "test.h5" 

         integer(hid_t) :: file_id, plist_id 
         integer :: error, mpierror

         call mpi_init(mpierror)
         call h5open_f(error) 

         call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
         call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

         call h5pclose_f(plist_id, error)
         call h5fclose_f(file_id, error)
         call h5close_f(error)

         call mpi_finalize(mpierror)

         end program])],
      [AC_MSG_RESULT([ok])], 
      [AC_MSG_FAILURE([Could not compile and link with parallel HDF5])])
  fi
  dnl Don't forget to pop
  AC_LANG_POP()

])
