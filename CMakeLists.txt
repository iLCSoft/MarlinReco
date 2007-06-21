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



####################################################################################################
# User section: modify as needed
####################################################################################################

# force cmake to use g77 fortran compiler
# FIXME: detect and test gfortran
SET( CMAKE_Fortran_COMPILER_INIT g77 )

# project name
PROJECT( MarlinReco C CXX Fortran )

# project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 3 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 1 )


### SETTINGS #################################################################

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( INSTALL_DOC "Set to ON to install Documentation" ON )

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
    SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" )
ENDIF()

### DOCUMENTATION ############################################################

FIND_PACKAGE( LATEX )
IF( LATEX_COMPILER AND MAKEINDEX_COMPILER AND PDFLATEX_COMPILER
    AND DVIPS_CONVERTER AND LATEX2HTML_CONVERTER )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${CMAKE_SOURCE_DIR}/doc/manual_html"
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${DVIPS_CONVERTER}" ARGS -q manual.dvi
        COMMAND "${LATEX2HTML_CONVERTER}" ARGS
            -mkdir -dir ../manual_html -verbosity 0 -info 0 -no_auto_link -split 0 -no_navigation manual.tex
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/doc/manual"
        COMMENT "Building Documentation - Users Manual..."
        VERBATIM )

    # add doc target
    ADD_CUSTOM_TARGET( doc DEPENDS
        "${CMAKE_SOURCE_DIR}/doc/manual_html" )
ELSE()
    MESSAGE( SEND_ERROR "Latex not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# install documentation
IF( INSTALL_DOC )
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${CMAKE_BINARY_DIR} ARGS doc)" )
    # install documentation
    INSTALL( DIRECTORY "${CMAKE_SOURCE_DIR}/doc"
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

# add library install path to the rpath list
SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH )

# append link pathes to rpath list
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH 1 )
MARK_AS_ADVANCED( CMAKE_INSTALL_RPATH_USE_LINK_PATH )

# output directories
SET( EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )

# DEPENDENCIES: this code has to be placed before adding any library or
# executable so that these are linked properly against the dependencies
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    # load macro
    IF( NOT EXISTS "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
        MESSAGE( FATAL_ERROR
            "\nSorry, could not find MacroCheckDeps.cmake...\n"
            "Please set CMAKE_MODULE_PATH correctly with: "
            "cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
    ENDIF()
    INCLUDE( "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
    CHECK_DEPS()
ENDIF()

# input directories
ADD_SUBDIRECTORY( src )

# create uninstall configuration file 
CONFIGURE_FILE( "${CMAKE_SOURCE_DIR}/cmake_uninstall.cmake.in"
                "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
                IMMEDIATE @ONLY )

# create uninstall target
ADD_CUSTOM_TARGET( uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake" )

# create configuration file from .in file
CONFIGURE_FILE( "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
                "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY )

# install configuration file
INSTALL( FILES "${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION . )

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
    "Set to ON to install Documentation" FORCE )

# export library dependencies (keep this as the last line in the file)
EXPORT_LIBRARY_DEPENDENCIES( "${PROJECT_NAME}LibDeps.cmake" )

