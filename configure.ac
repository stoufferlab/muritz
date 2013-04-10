
AC_INIT([muritz], [0.01], [daniel.stouffer@canterbury.ac.nz])
AC_PREREQ(2.61)
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#
# Init Automake subsystem
#
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])


# Indicate this is a release build and that
# dependancies for changes between Makefile.am and Makefile.in
# should not be checked. This makes compiling faster. If you are working
# on the library, run: ./configure --enable-maintainer-mode to enable the dependancies
#dnl AM_MAINTAINER_MODE

# Need the C++ Compiler, Installer and ranlib to strip libraries
#
#AC_PROG_CXX
#AC_LANG(C++)
#AC_PROG_INSTALL
#AC_PROG_LN_S
#AC_TYPE_SIZE_T
#AC_HEADER_STDBOOL
#AC_FUNC_STRTOD
#AC_PROG_CC
#AC_PROG_LIBTOOL

# Check for GSL
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
