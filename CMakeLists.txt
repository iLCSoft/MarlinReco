########################################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6.2 FATAL_ERROR )
########################################################



# due to the incompatibility between g77/gfortran on 32/64bit machines
# it is enforced to set the fortran compiler before running cmake
IF( NOT DEFINED ENV{FC} )
    MESSAGE( FATAL_ERROR "please set FC before running cmake, e.g.: export FC=g77" )
ENDIF()

# fix for lcio: all c++ sources need to be compiled with -Wno-long-long
SET( ENV{CXXFLAGS} "-Wno-long-long $ENV{CXXFLAGS}" )

# require proper c and c++
# -pedantic cannot be set globally because it cannot be removed for certain
# source files only: grep pedantic src/CMakeLists.txt
SET( ENV{CXXFLAGS} "-Wall -ansi $ENV{CXXFLAGS}" ) # FIXME -pedantic
SET( ENV{CFLAGS} "-Wall -ansi -pedantic $ENV{CFLAGS}" )

####################################################################################################
# User section: modify as needed
####################################################################################################


# project name
PROJECT( MarlinReco C CXX Fortran )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 19 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### SETTINGS ################################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )




SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DNDEBUG" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -DNDEBUG" )

GET_FILENAME_COMPONENT( FORTRAN_COMPILER "${CMAKE_Fortran_COMPILER}" NAME_WE )
IF( FORTRAN_COMPILER STREQUAL "g77" )
    SET( FORTRAN_LIBNAME "g2c" )
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-globals" ) # -ff2c ?
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dg77Fortran" ) # -Df2cFortran

ELSE() # gfortran / gfortran-32

    SET( FORTRAN_LIBNAME "gfortran" )
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DgFortran" )
ENDIF()

IF( NOT BUILD_32BIT_COMPATIBLE )
    SET( FORTRAN_LIB_PATHS
        # 64 bit paths
        /usr/lib/gcc/x86_64-redhat-linux/3.4.3      # SL4 native 64bit
        /usr/lib/gcc/x86_64-redhat-linux/3.4.6      # SL5 native 64bit
        /usr/lib/gcc/x86_64-redhat-linux/4.1.1      # SL5 native 64bit
        /usr/lib64
    )
ENDIF()

SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        # 32 bit paths
        /usr/lib/gcc-lib/i386-redhat-linux/3.2.3    # SL3 native 32bit
        /usr/lib/gcc/i386-redhat-linux/3.4.3        # SL4 native 32bit
        /usr/lib/gcc/i386-redhat-linux/3.4.6        # SL5 native 32bit
        /usr/lib/gcc/i386-redhat-linux/4.1.1        # SL5 native 32bit
)

IF( BUILD_32BIT_COMPATIBLE )
    SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        /usr/lib/gcc/x86_64-redhat-linux/3.4.3/32   # SL4 64bit 32bit compatible
        /usr/lib/gcc/x86_64-redhat-linux/3.4.6/32   # SL5 64bit 32bit compatible
        /usr/lib/gcc/x86_64-redhat-linux/4.1.1/32   # SL5 64bit 32bit compatible
        /usr/lib
    )
ENDIF()

SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        /usr/lib/gcc/x86_64-linux-gnu/4.4           # Ubuntu 10.04 64bit
        /usr/lib/gcc/i486-linux-gnu/4.4             # Ubuntu 10.04 32bit
        /usr/osxws/lib                              # OSX
)

#MESSAGE( STATUS "FORTRAN_LIB_PATHS: ${FORTRAN_LIB_PATHS}" )

SET( FORTRAN_LIB FORTRAN_LIB-NOTFOUND )
FIND_LIBRARY( FORTRAN_LIB NAMES ${FORTRAN_LIBNAME} PATHS ${FORTRAN_LIB_PATHS} NO_DEFAULT_PATH )
# if previous find fails, try to find in system paths..
FIND_LIBRARY( FORTRAN_LIB NAMES ${FORTRAN_LIBNAME} )

IF( NOT FORTRAN_LIB )
    MESSAGE( FATAL_ERROR "failed to find ${FORTRAN_LIBNAME} library!!" )
ELSE()
    MESSAGE( STATUS "Check for lib${FORTRAN_LIBNAME}: ${FORTRAN_LIB}" )
ENDIF()
MARK_AS_ADVANCED( FORTRAN_LIB )



### DEPENDENCIES ############################################################


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( GSL REQUIRED )
FIND_PACKAGE( CERNLIB REQUIRED )
FIND_PACKAGE( RAIDA REQUIRED )
FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathMore )

FOREACH( pkg Marlin MarlinUtil GSL CERNLIB RAIDA ROOT )
    # fix for backwards compatibility
    IF( NOT ${pkg}_DIR AND DEFINED ${pkg}_HOME )
        SET( ${pkg}_DIR ${${pkg}_HOME} )
    ENDIF()

    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

# MathMore library from ROOT is required for BCalTagEfficiency
LINK_LIBRARIES( ${ROOT_MATHMORE_LIBRARY} )
# ============================================================================



### DOCUMENTATION ############################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )


ADD_CUSTOM_TARGET( doc )

# MANUAL
FIND_PACKAGE( LATEX )
IF( LATEX_COMPILER AND MAKEINDEX_COMPILER AND PDFLATEX_COMPILER
    AND DVIPS_CONVERTER AND LATEX2HTML_CONVERTER )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/manual_html"
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${DVIPS_CONVERTER}" ARGS -q -o manual.ps manual.dvi
        COMMAND "${LATEX2HTML_CONVERTER}" ARGS
            -mkdir -dir ../manual_html -verbosity 0 -info 0 -no_auto_link -split 0 -no_navigation manual.tex
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc/manual"
        COMMENT "Building Documentation - Users Manual..."
        VERBATIM )

    # add doc target
    ADD_CUSTOM_TARGET( manual DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/manual_html" )
    # tell doc target to call manual
    ADD_DEPENDENCIES( doc manual )
ELSE()
    MESSAGE( STATUS "Latex not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()


# APIDOC
FIND_PACKAGE( Doxygen )
IF( DOXYGEN_FOUND )
                                                                                                                                                            
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/html"
#        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
#            "${PROJECT_SOURCE_DIR}/README" "${PROJECT_SOURCE_DIR}/doc/README"
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS "Doxyfile"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation..."
        VERBATIM )
                                                                                                                                                            
    # add doc target
    ADD_CUSTOM_TARGET( apidoc DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/html" )
    # tell doc target to call apidoc
    ADD_DEPENDENCIES( doc apidoc )
ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()
    

# install documentation
IF( INSTALL_DOC )

    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_BUILD_TOOL} doc)" )

    # install documentation
    INSTALL_DIRECTORY( "${PROJECT_SOURCE_DIR}/doc" DESTINATION . )
ENDIF()


# input directories
ADD_SUBDIRECTORY( src )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )

