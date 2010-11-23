########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
########################################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.6.2 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
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
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 19 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


### SETTINGS #################################################################


# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( BUILD_32BIT_COMPATIBLE "Set to ON to build in 32 bit compatibility mode" OFF )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )


SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DNDEBUG" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -DNDEBUG" )

GET_FILENAME_COMPONENT( FORTRAN_COMPILER "${CMAKE_Fortran_COMPILER}" NAME_WE )
IF( FORTRAN_COMPILER STREQUAL "g77" )
    SET( FORTRAN_LIBNAME "g2c" )
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-globals" ) # -ff2c ?
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dg77Fortran" ) # -Df2cFortran

ENDIF()
IF( FORTRAN_COMPILER STREQUAL "gfortran" )
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


# ========= project dependencies =============================================
# MathMore library from ROOT is required for BCalTagEfficiency
SET( ROOT_USE_COMPONENTS MathMore )
SET( req_pkgs Marlin MarlinUtil GEAR LCIO CLHEP GSL CERNLIB RAIDA ROOT )
IF(APPLE)
    LIST( APPEND req_pkgs CED )
ENDIF(APPLE)
FOREACH( pkg ${req_pkgs} )
    # fix for backwards compatibility
    IF( NOT ${pkg}_DIR AND DEFINED ${pkg}_HOME )
        SET( ${pkg}_DIR ${${pkg}_HOME} )
    ENDIF()

    FIND_PACKAGE( ${pkg} REQUIRED )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
ENDFOREACH()

LINK_LIBRARIES( ${ROOT_MATHMORE_LIBRARY} )
# ============================================================================

# set default cmake build type to RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()


#---------------- 32/64 bit issues ---------------------------------------
IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "64 bit architecture detected" )

    IF( BUILD_32BIT_COMPATIBLE )
        IF( COMMAND SET_PROPERTY )
            SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS 0)
            SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32" )
            SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32" )
            SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m32" )
            SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32" )
            SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -m32" )
            SET( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -m32" )
        ELSE()
            MESSAGE( FATAL_ERROR "to build in 32 bit compatibility mode you need cmake >= 2.6" )
        ENDIF()
        MESSAGE( STATUS "Build in 32 bit compatibility mode" )
    ELSE()
        MESSAGE( STATUS "Build in native 64 bit mode" )
    ENDIF()
ENDIF()
#-------------------------------------------------------------------------


### DOCUMENTATION ############################################################

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
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${PROJECT_BINARY_DIR} ARGS doc)" )
    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE
            PATTERN "*/.svn*" EXCLUDE
            )
ENDIF()


##########################################################################################
# End of User section: please try not to modify below this line
##########################################################################################


# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

# add install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# add install path to the rpath list (apple)
IF( APPLE )
    SET( CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" )
    MARK_AS_ADVANCED( CMAKE_INSTALL_NAME_DIR )
ENDIF()

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# input directories
ADD_SUBDIRECTORY( src )

# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION . )

# display status message for important variables
MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}" )
#IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
#    MESSAGE( STATUS "BUILD_32BIT_COMPATIBLE = ${BUILD_32BIT_COMPATIBLE}" )
#ENDIF()
MESSAGE( STATUS "INSTALL_DOC = ${INSTALL_DOC}" )
MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

# force some variables that could be defined in the command line
# to be written to cache
SET( BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL
    "Set to OFF to build static libraries" FORCE )
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
    "Where to install ${PROJECT_NAME}" FORCE )
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
SET( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" CACHE PATH
    "Path to custom CMake Modules" FORCE )
SET( INSTALL_DOC "${INSTALL_DOC}" CACHE BOOL
    "Set to OFF to skip build/install Documentation" FORCE )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )
INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}LibDeps.cmake" DESTINATION lib/cmake )
