########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
########################################################


########################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
#SET( CMAKE_BACKWARDS_COMPATIBILITY 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
# allow more human readable "if then else" constructs
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
########################################################


# set fortran compiler
SET( CMAKE_Fortran_COMPILER_INIT g77 )

####################################################################################################
# User section: modify as needed
####################################################################################################

# project name
PROJECT( MarlinReco C CXX Fortran )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 10 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 4 )


### SETTINGS #################################################################

# fortran flags
SET( CMAKE_Fortran_FLAGS "-ff2c" CACHE STRING
    "Flags for Fortran compiler." )
# gfortran does not support -Wno-globals
GET_FILENAME_COMPONENT( FORTRAN_COMPILER "${CMAKE_Fortran_COMPILER}" NAME_WE )
IF( FORTRAN_COMPILER STREQUAL "g77" )
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-globals" )
ENDIF()
SET( CMAKE_Fortran_FLAGS_DEBUG "-g" CACHE STRING
    "Flags used by the compiler during debug builds." )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "-Os -DNDEBUG" CACHE STRING
    "Flags used by the compiler during release minsize builds." )
SET( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING
    "Flags used by the compiler during Release with Debug Info builds." )
SET( CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files" )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

# project dependencies
# e.g. SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil LCIO GEAR CLHEP GSL RAIDA" )
SET( ${PROJECT_NAME}_DEPENDS "Marlin MarlinUtil GEAR LCIO CLHEP GSL CERNLIB" )

# set default cmake build type to RelWithDebInfo
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE "RelWithDebInfo" )
ENDIF()

# set default install prefix to project root directory
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

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
            PATTERN "*CVS*" EXCLUDE )
ENDIF()


##########################################################################################
# End of User section: please try not to modify below this line
##########################################################################################


# library *nix style versioning
SET( ${PROJECT_NAME}_VERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )
SET( ${PROJECT_NAME}_SOVERSION
    "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}" )

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

# DEPENDENCIES: this code has to be placed before adding any library or
# executable so that these are linked properly against the dependencies
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    IF( NOT DEFINED CMAKE_MODULE_PATH )
        MESSAGE( FATAL_ERROR "CMAKE_MODULE_PATH not set! Set it with: "
            "-DCMAKE_MODULE_PATH=\"/path_to_ilcsoft/CMakeModules\"" )
    ENDIF()
    # load macro
    INCLUDE( "MacroCheckDeps" )
    CHECK_DEPS()
ENDIF()

# input directories
ADD_SUBDIRECTORY( src )

# create uninstall configuration file 
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall
  "${CMAKE_COMMAND}" -P "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake" )

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
MESSAGE( STATUS "${PROJECT_NAME}_DEPENDS = \"${${PROJECT_NAME}_DEPENDS}\"" )
MESSAGE( STATUS "BUILD_WITH = \"${BUILD_WITH}\"" )
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

# export build settings
INCLUDE( CMakeExportBuildSettings )
CMAKE_EXPORT_BUILD_SETTINGS( "${PROJECT_NAME}BuildSettings.cmake" )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )

