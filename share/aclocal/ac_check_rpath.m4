dnl Checking the compiler option which is used to set the
dnl runtime library path.

AC_DEFUN([AC_CHECK_RPATH_OPTION], [
AC_MSG_CHECKING(whether to use '-R' or '-Wl,-rpath')
AC_CACHE_VAL(ac_cv_check_rpath_option,
[
testfile='.ac_cv_check_rpath_option_test_file.c'
echo 'int main() {}' > $testfile
for i in '-R' '-Wl,-rpath,'; do
    if $CC ${i}/lib -o /dev/null $testfile 2> /dev/null;
    then ac_cv_check_rpath_option=$i;
    fi
done
\rm -rf $testfile
])
if test -z "$ac_cv_check_rpath_option"; then
    AC_MSG_WARN(Could not find a compiler option to specify the run time library directory (like '-R' or '-Wl,-rpath ') for the compiler '$CC'. In consequence the variable 'RPATH' will remain unset.)
    RPATH=@RPATH@
else
    RPATH=$ac_cv_check_rpath_option
    AC_MSG_RESULT($RPATH)
fi
AC_SUBST(RPATH)
])
