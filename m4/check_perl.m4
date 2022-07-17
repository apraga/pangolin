# AC_CHECK_PERLMODULE([module])
# Search if the module is installed
AC_DEFUN([AC_CHECK_PERLMODULE],
[
  AC_MSG_CHECKING(for Perl module $1)
  # Must redirect stderr
  out=`perldoc -l $1 2>&1`
  # Must use silent option for grep
  if echo $out | grep -q "No documentation" ; then
    AC_MSG_RESULT([failed, will not be able to run test suite])
  else
    AC_MSG_RESULT(ok)
  fi
])
