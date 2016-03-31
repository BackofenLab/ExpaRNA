#
AC_INIT([ExpaRNA], [1.0], [heyne@informatik.uni-freiburg.de])
AC_CONFIG_AUX_DIR(config)

# store user given flags for possible later resetting
TMP_CXXFLAGS=$CXXFLAGS
TMP_CPPFLAGS=$CPPFLAGS
TMP_CFLAGS=$CFLAGS

# compiler check (if CXXFLAGS or CFLAGS are not available '-g -O2' is set)
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# enable/disable debugging information if wanted
AC_ARG_ENABLE(debug-info,
	AC_HELP_STRING(
	    [--enable-debug-info],
	    [enable '-g' compiler flag for debug information build (default=disabled)]
	),
	debuginfo="$enableval",
  	debuginfo=no
)

# test if debugging information is wanted and handle setting
AC_MSG_CHECKING(whether to build with debug information (bigger executable))
AC_MSG_RESULT($debuginfo)
if test  x"$debuginfo" = x"no"; then
	# undo '-g -O2' debuggin add of AC_PROG_CXX and AC_PROG_CC macro, if flags are empty
	#CXXFLAGS="$TMP_CXXFLAGS -std=gnu++0x"
	#CPPFLAGS="$TMP_CPPFLAGS -std=gnu++0x"
	#CFLAGS="$TMP_CFLAGS -std=gnu++0x"
	CXXFLAGS=$TMP_CXXFLAGS
	CPPFLAGS=$TMP_CPPFLAGS
	CFLAGS=$TMP_CFLAGS
else
	CXXFLAGS="$CXXFLAGS -g"
	CPPFLAGS="$CPPFLAGS -g"
	CFLAGS="$CFLAGS -g"
fi


AC_PROG_RANLIB
AM_INIT_AUTOMAKE
AC_LANG_CPLUSPLUS

###############################################################################
# BEGIN STATIC LINKING CHECK
###############################################################################

AC_ARG_ENABLE(static-linking,
    AC_HELP_STRING(
    	[--enable-static-linking],
	    [enables static linking of the executables (default=disabled)]
	),
    [enable_static_linking=$enableval],
    [enable_static_linking=no]
)

if test "$enable_static_linking" = "yes"; then
	LDFLAGS="-static $LDFLAGS"
#	LIBS="-lpthread $LIBS"
fi
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT($enable_static_linking)

###############################################################################
# END STATIC LINKING CHECK
###############################################################################


###############################################################################
# BEGIN PATH-SUPPORT CHECK
###############################################################################

# Vienna RNA package library path support, if not installed in usual directories
LIBS="-lRNA $LIBS"
AC_ARG_WITH(RNA,
    AC_HELP_STRING(
        [--with-RNA=prefix],
        [alternative prefix path to Vienna RNA library]
      ),
      RNAPATHSET=1,
      RNAPATHSET=0
)
if test  $RNAPATHSET = 1 ; then
        CPPFLAGS="$CPPFLAGS -I$with_RNA/include"
        CFLAGS="$CFLAGS -I$with_RNA/include"
        LDFLAGS="-L$with_RNA/lib $LDFLAGS"
fi

# xerces-c Library path support, if not installed in usual directories
# LIBS="-lxerces-c $LIBS"
# AC_ARG_WITH(xml,
#    AC_HELP_STRING(
#        [--with-xml=prefix],
#        [alternative prefix path to xerces-c library]
#      ),
#      xmlPATHSET=1,
#      xmlPATHSET=0
#)		
#if test  $xmlPATHSET = 1 ; then
#        CPPFLAGS="$CPPFLAGS -I$with_xml/include"
#        CFLAGS="$CFLAGS -I$with_xml/include"
#        LDFLAGS="-L$with_xml/lib $LDFLAGS"
#fi

# LIBS="-lleda $LIBS"

###############################################################################
# END PATH-SUPPORT CHECK
###############################################################################

LDFLAGS="$LDFLAGS $LIBS"
AC_SUBST(GLOBALCPPFLAGS, $CPPFLAGS )
AC_SUBST(GLOBALLDFLAGS, $LDFLAGS )

#
# Checks for libraries.
#
AC_CHECK_LIB(RNA, [energy_of_struct], [xmlNOTFOUND=0], [ xmlNOTFOUND=1])

###############################################################################
# BEGIN VIENNA CHECK
###############################################################################


# check for Vienna RNA headers
AC_MSG_CHECKING([for the Vienna RNA C header])
AC_COMPILE_IFELSE(
     [ #include <ViennaRNA/fold_vars.h> ],
     [
        AC_MSG_RESULT(yes)
        AC_CHECK_LIB(RNA, [fold], [RNANOTFOUND=0], [ RNANOTFOUND=1])
         ],
     [
       AC_MSG_RESULT(no)
       RNANOTFOUND=1;
     ]
)

###############################################################################
# END VIENNA CHECK
###############################################################################

###############################################################################
# BEGIN xerces-c CHECK
###############################################################################


# check for xerces-c headers
#AC_MSG_CHECKING([for the xerces-c header])
#AC_COMPILE_IFELSE(
#     [ #include <xercesc/util/XercesVersion.hpp> ],
#     [	
#        AC_MSG_RESULT(yes)
#	xmlNOTFOUND=0;
#        #AC_CHECK_LIB(xerces-c, [XMLPlatformUtils::Initialize], [xmlNOTFOUND=0], [ xmlNOTFOUND=1])
#         ],
#     [
#       AC_MSG_RESULT(no)
#       xmlNOTFOUND=1;
#     ]
#)

###############################################################################
# END xerces-c CHECK
###############################################################################
#AX_LIB_XERCES([2.7.0])

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
# BEGIN DEBUG/PROFILING CHECK
###############################################################################
AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (default=disabled)]
	),
 	debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE(DEBUG)
	CPPFLAGS="$CPPFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
	LDFLAGS="$LDFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
else
	AC_DEFINE(NDEBUG)
	CPPFLAGS="$CPPFLAGS -O2"
	# add -std=gnu++0x for new standard, unorderd_map instead of ext/hash_map
fi

AC_MSG_CHECKING(whether to build with profiler information)
profiler=no
AC_ARG_ENABLE(profiling,
	AC_HELP_STRING(
	    [--enable-profiling],
	    [enable profiling data generation (default=disabled)]
	),
 	profiler="$enableval"
)
AC_MSG_RESULT($profiler)
if test x"$profiler" = x"yes"; then
	CPPFLAGS="$CPPFLAGS -pg"
    LDFLAGS="$LDFLAGS -pg"
fi

###############################################################################
# END DEBUG/PROFILING CHECK
###############################################################################


###############################################################################
# BEGIN DO NOT COMPILE CHECKS
###############################################################################

DEPENDENCYNOTFOUND=0;

# error output if ViennaRNA not found
if test "$RNANOTFOUND" = "1"; then
    AC_MSG_NOTICE()
    AC_MSG_NOTICE(the Vienna RNA C library is required for building.)
    AC_MSG_NOTICE(it can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
    AC_MSG_NOTICE()
     if test "$RNAPATHSET" = "1"; then
        AC_MSG_NOTICE(cant find Vienna RNA library in given path '$with_RNA'.)
     else
        AC_MSG_NOTICE(if Vienna RNA is installed elsewhere use --with-RNA=prefix.)
     fi
     DEPENDENCYNOTFOUND=1;
fi


# error output if xerces-c not found
#if test "$HAVE_XERCES" = "0"; then
#    AC_MSG_NOTICE()
#    AC_MSG_NOTICE(the xerces-c library is required for building.)
#    AC_MSG_NOTICE(cant find xerces-c library in given path '$with_xerces'.)
#    DEPENDENCYNOTFOUND=1;
#fi


# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
    AC_MSG_NOTICE()
    AC_MSG_ERROR(library not found! Try './configure --help'.)
fi

###############################################################################
# END DO NOT COMPILE CHECKS
###############################################################################


AC_CONFIG_FILES([ Makefile  src/Makefile ])
#
#
AC_OUTPUT
