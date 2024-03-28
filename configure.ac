#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([TSEM], [1.1.0], [])
AM_INIT_AUTOMAKE([subdir-objects tar-ustar])
LT_INIT
AC_CONFIG_SRCDIR([examples/poisson_mixBC/poisson.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX([mpicxx])
AC_PROG_CC([mpicc])
AC_PROG_LN_S
dnl  AC_PROG_RANLIB

AC_LANG(C++)

AC_ARG_WITH([afepack], [AS_HELP_STRING([--with-afepack=PATH], [Specify the installation path for AFEPACK])], [AFEPackPath="$withval"], [AFEPackPath="/usr/local"])
AC_SUBST([AFEPackPath])

AC_ARG_WITH([dealii], [AS_HELP_STRING([--with-dealii=PATH], [Specify the installation path for deal.II])], [dealIIPath="$withval"], [dealIIPath="/usr/local"])
AC_SUBST([dealIIPath])

# Checks for libraries.
# FIXME: Replace 'main' with a function in '-lAFEPack':
AC_CHECK_LIB([AFEPack], [main])
# FIXME: Replace 'main' with a function in '-ldeal_II':
AC_CHECK_LIB([deal_II], [main],
  [deal_II_library=yes, LIBS="-ldeal_II -L\$(AFEPackPath)/lib -L\$(dealIIPath)/lib $LIBS", CPPFLAGS="-I\$(AFEPackPath)/include -I\$(dealIIPath)/include/deal.II -I\$(dealIIPath)/include $CPPFLAGS"],
  [deal_II_library=no])
 if test "x$deal_II_library" == xno; then
  AC_CHECK_LIB([base], [main],
    [deal_II_library=yes, LIBS="-llac -lbase $LIBS"],
    [deal_II_library=no])
  if test "x$deal_II_library" == xno; then
    AC_MSG_ERROR(deal.II library not found)
  fi
 fi

# FIXME: Replace 'main' with a function in '-ltbb':
AC_CHECK_LIB([tbb], [main],
  [tbb_library=yes, LIBS="-ltbb $LIBS"],
  [tbb_library=no])



# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_RESTRICT
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_SUBST([AM_CPPFLAGS],['-I$(top_srcdir)/include'])
AC_SUBST([AM_CXXFLAGS],['-D__SERIALIZATION -DMULTITHREAD -pthread -fno-delete-null-pointer-checks -fpermissive -llapack -llapacke'])

AC_CONFIG_FILES([Makefile
		 include/Makefile
		 examples/Makefile
		 examples/KohnSham-Helium/Makefile
		 examples/laplaceEig/Makefile
		 examples/movingMesh/Makefile
		 examples/poisson_mixBC/Makefile
                 src/Makefile])
AC_OUTPUT
