#############################################################################
# EXAMPLE cmake build setup for MarlinReco
#
# For building MarlinUtil with cmake:
#
# (0) --  modify the pathes and versions of the software packages !!!!
#
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
# @author Jan Engels, DESY
#############################################################################

#############################################################################
# Setup path variables  <<<< MODIFY ACCORDINGLY !!!!!!
#############################################################################

# ILC_HOME
SET( ILC_HOME "/afs/desy.de/group/it/ilcsoft/v01-03"
    CACHE PATH "Path to ILC Software" FORCE )

# Path to Marlin
SET( Marlin_HOME "${ILC_HOME}/Marlin/v00-09-10"
    CACHE PATH "Path to Marlin" FORCE )

# Path to MarlinUtil
SET( MarlinUtil_HOME "${ILC_HOME}/MarlinUtil/v00-05"
    CACHE PATH "Path to MarlinUtil" FORCE )

# Path to LCIO
SET( LCIO_HOME "${ILC_HOME}/lcio/v01-09"
    CACHE PATH "Path to LCIO" FORCE )

# Path to GEAR
SET( GEAR_HOME "${ILC_HOME}/gear/v00-08"
    CACHE PATH "Path to GEAR" FORCE )

# Path to CLHEP
SET( CLHEP_HOME "${ILC_HOME}/CLHEP/2.0.3.1"
    CACHE PATH "Path to CLHEP" FORCE )

# Path to GSL
SET( GSL_HOME "${ILC_HOME}/gsl/1.8"
    CACHE PATH "Path to GSL" FORCE )

# Path to CERNLIB
SET( CERNLIB_HOME "${ILC_HOME}/cernlib/2006"
    CACHE PATH "Path to CERNLIB" FORCE )

# CMake Modules Path
SET( CMAKE_MODULE_PATH "${ILC_HOME}/CMakeModules/v01-05"
    CACHE PATH "Path to CMake Modules" FORCE )

###############################################
# Project Options
###############################################

#SET( INSTALL_DOC OFF CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )

# set cmake build type
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
#SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE )

###############################################
# Advanced Options
###############################################

#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

# installation path for MarlinReco
#SET( CMAKE_INSTALL_PREFIX "/foo/bar" CACHE STRING "Where to install MarlinReco" FORCE )
