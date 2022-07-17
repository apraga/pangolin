# GET_MPI_COMPILER([possible names], [compiler])
# Search in PATH all mpi wrapper and select the correct one
# The configure script already defines PATH_SEPARATOR, which allow us to parse
# PATH
# Tests all the executable defined
AC_DEFUN([GET_MPI_COMPILER],
[
  as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
  for as_dir in $PATH
  do
     IFS=$as_save_IFS
     test -z "$as_dir" && as_dir=.
     for ac_exec_ext in '' $1; do
     if test -f "$as_dir/$ac_exec_ext"; then
        GET_MPI_COMPILER_VERSION($as_dir"/"$ac_exec_ext, mpi_version)
        if test x"$mpi_version" = x"$2"; then
           MPI_DIR=$as_dir
           MPI_WRAPPER=$ac_exec_ext
           break 2
        fi
     fi
     done
  done
  IFS=$as_save_IFS

  # If not found, use default
  if test x"$MPI_DIR" != x; then
  # Remove trailing slashes
  # Use quadrigraph @<:@ for escaping the [ character
    MPI_DIR=`expr "X$MPI_DIR" : 'X\(.*@<:@^/]\)' \| "X$MPI_DIR" : 'X\(.*\)'`
  fi
])

# Get the result of mpif90 -show (name of the compiler) to mpi_version
# GET_MPI_COMPILER_VERSION([exec name, output])
AC_DEFUN([GET_MPI_COMPILER_VERSION],
[
    $2=`$1 -show | awk '{ print $ 1}'`
])
