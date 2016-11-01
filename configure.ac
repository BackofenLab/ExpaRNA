
AC_PREREQ([2.65])
# 5 argument version only available with aclocal >= 2.64
AC_INIT( [ExpaRNA], [1.0.1], [], [expaRNA], [http://www.bioinf.uni-freiburg.de] )

AC_CONFIG_AUX_DIR([.])
AC_CONFIG_SRCDIR([src/Util_String.cpp])

# store user given flags for possible later resetting
TMP_CXXFLAGS=$CXXFLAGS
TMP_CPPFLAGS=$CPPFLAGS
TMP_CFLAGS=$CFLAGS

# compiler check (if CXXFLAGS or CFLAGS are not available '-g -O2' is set)
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB


# automake initialisation (mandatory) and check for minimal automake API version
AM_INIT_AUTOMAKE([1.11])

# use the C++ compiler for the following checks
AC_LANG([C++])


DEPENDENCYNOTFOUND=0;


###############################################################################
# DEBUG SUPPORT SETUP
###############################################################################

AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	    [enable debug data generation (def=no)])],
  	[debuger="$enableval"])
AC_MSG_RESULT([$debuger])
if test x"$debuger" = x"yes"; then
	AC_DEFINE([_DEBUG], [1], [Run in DEBUG mode with additional assertions and debug output])
	AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
	AM_CFLAGS="$AM_CFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
else
	AC_DEFINE([NDEBUG], [1], [Run in normal mode with minimal assertions])
	AM_CXXFLAGS="$AM_CXXFLAGS -O3 -Wno-deprecated"
	AM_CFLAGS="$AM_CFLAGS -O3 -Wno-deprecated"
fi


###############################################################################
# BEGIN STATIC LINKING CHECK
###############################################################################

AC_ARG_ENABLE([static-linking],
    [AC_HELP_STRING(
    	[--enable-static-linking],
	    [enables static linking of the executables (default=disabled)]
	)],
    [enable_static_linking=$enableval],
    [enable_static_linking=no]
)
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT([$enable_static_linking])

if test "$enable_static_linking" = "yes"; then
	AM_LDFLAGS="-static $AM_LDFLAGS"
#	LIBS="-lpthread $LIBS"
fi

###############################################################################
# END STATIC LINKING CHECK
###############################################################################


###############################################################################
# BEGIN PROFILING CHECK
###############################################################################

AC_MSG_CHECKING([whether to build with profiler information])
profiler=no
AC_ARG_ENABLE([profiling],
	[AC_HELP_STRING(
	    [--enable-profiling],
	    [enable profiling data generation (default=disabled)]
	)],
 	[profiler="$enableval"]
)
AC_MSG_RESULT([$profiler])
if test x"$profiler" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -pg"
    AM_LDFLAGS="$AM_LDFLAGS -pg"
fi

###############################################################################
# Vienna RNA package library path support, if not installed in usual directories
###############################################################################
AC_ARG_WITH([RNA],
	[AC_HELP_STRING(
	    [--with-RNA=PREFIX],
	    [alternative prefix path to Vienna RNA library]
	  )],
	  [RNAPATHSET=1],
	  [RNAPATHSET=0]
)
if test  $RNAPATHSET = 1 ; then
	# set compiler and linker flags if needed
	AM_CXXFLAGS="-I$with_RNA/include $AM_CXXFLAGS -fopenmp"
	AM_CFLAGS="-I$with_RNA/include $AM_CFLAGS -fopenmp"
	AM_LDFLAGS="-L$with_RNA/lib $AM_LDFLAGS"
fi

###############################################################################
# BEGIN VIENNA CHECK
###############################################################################
# check for Vienna RNA headers
AC_MSG_CHECKING([for the Vienna RNA package headers])
OLD_CPPFLAGS=$CPPFLAGS
OLD_CXXFLAGS=$CXXFLAGS
OLD_LDFLAGS=$LDFLAGS
CPPFLAGS="$CPPFLAGS $AM_CXXFLAGS"
CXXFLAGS="$CXXFLAGS $AM_CXXFLAGS"
LDFLAGS="$LDFLAGS $AM_LDFLAGS"
AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM([[#include <ViennaRNA/fold_vars.h>]],[[int var = noGU;]])],
     [ 
		AC_MSG_RESULT([yes])
#		AC_CHECK_LIB([RNA], [var_time_stamp], [RNANOTFOUND=0], [ RNANOTFOUND=1])
     ],
     [
       AC_MSG_RESULT([no])
	   AC_MSG_NOTICE([DEBUG : used CXXFLAGS = $CXXFLAGS])
	   AC_MSG_NOTICE([DEBUG : used CPPFLAGS = $CPPFLAGS])
	   AC_MSG_NOTICE([DEBUG : used  LDFLAGS = $LDFLAGS])
       RNANOTFOUND=1;
     ]
)
CPPFLAGS=$OLD_CPPFLAGS
CXXFLAGS=$OLD_CXXFLAGS
LDFLAGS=$OLD_LDFLAGS

# error output if ViennaRNA not found
if test "$RNANOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE([The Vienna RNA C library is required.])
	AC_MSG_NOTICE([ -> It can be obtained from http://www.tbi.univie.ac.at/.])
	AC_MSG_NOTICE()
 	if test "$RNAPATHSET" = "1"; then
		AC_MSG_NOTICE([ -> Can't find the Vienna RNA library in given path '$with_RNA'.])
 	else
		AC_MSG_NOTICE([ -> If installed in a non-standard path, please use '--with-RNA=PREFIX'.])
 	fi
 	DEPENDENCYNOTFOUND=1;
else
	# register Vienna RNA lib for linking
	LIBS="$LIBS -lRNA"
fi


###############################################################################
# END VIENNA CHECK
###############################################################################



###############################################################################
# END PATH-SUPPORT CHECK
###############################################################################

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string fstream ext/algorithm sys/time.h ext/hash_map vector map deque])
#
# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T
#

###############################################################################
# FINAL DEPENDENCY CHECK AND EXIT IF NEEDED
###############################################################################

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_ERROR([Some dependency was not met! See above for errors and relate to './configure --help'.])
fi

##########################################################################

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# setup files to create
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

# create files
AC_OUTPUT
