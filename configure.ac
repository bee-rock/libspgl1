AC_INIT([spgl1++],[0.0.2],[brock.hargreaves@gmail.com],[spgl1++],[http://sourceforge.net/projects/libspgl1/])
AC_CONFIG_SRCDIR([src/lib/spgl1++/all.hpp])
AC_CONFIG_HEADERS([src/lib/config.hpp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.9])
LT_INIT([win32-dll])
: ${CXXFLAGS=""}
AC_PROG_CXX
AC_LANG([C++])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
ACX_PTHREAD
LIBS="$PTHREAD_LIBS $LIBS"
AM_CXXFLAGS="$AM_CXXFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"

# debug compilation support
AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug build [default=no]])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Wextra -Weffc++ -pedantic"
else
    AM_CXXFLAGS="$AM_CXXFLAGS -O2"
fi

# distribute additional compiler and linker flags among Makefiles
# --> set and change these variables instead of CXXFLAGS or LDFLAGS (for user only)
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/lib/Makefile])
AC_CONFIG_FILES([test/Makefile])

AC_PROG_INSTALL

AC_OUTPUT
