########################################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6.2 FATAL_ERROR )
########################################################



# due to the incompatibility between g77/gfortran on 32/64bit machines
# it is enforced to set the fortran compiler before running cmake
#IF( NOT DEFINED ENV{FC} )
#    MESSAGE( FATAL_ERROR "please set FC before running cmake, e.g.: export FC=gfortran" )
#ENDIF()

# fix for lcio: all c++ sources need to be compiled with -Wno-long-long
SET( ENV{CXXFLAGS} "-Wno-long-long $ENV{CXXFLAGS}" )

# require proper c and c++
# -pedantic cannot be set globally because it cannot be removed for certain
# source files only: grep pedantic src/CMakeLists.txt
SET( ENV{CXXFLAGS} "-Wall -ansi $ENV{CXXFLAGS}" ) # FIXME -pedantic
SET( ENV{CFLAGS} "-Wall -ansi -pedantic $ENV{CFLAGS}" )


####################################################################################################


# project name
PROJECT( MarlinReco C CXX Fortran )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 19 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( GSL REQUIRED )
FIND_PACKAGE( CERNLIB REQUIRED )
FIND_PACKAGE( ROOT 5.27 REQUIRED COMPONENTS MathMore )

OPTION( MARLINRECO_AIDA "Set to ON to build MarlinReco with AIDA" ON )
IF( MARLINRECO_AIDA )
    FIND_PACKAGE( AIDA )
    IF( NOT AIDA_FOUND )
        MESSAGE( SEND_ERROR "AIDA not found. Please set MARLINRECO_AIDA to OFF or AIDA_DIR=/path/to/AIDA" )
    ENDIF()
ENDIF()

FOREACH( pkg Marlin MarlinUtil GSL CERNLIB AIDA ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

# MathMore library from ROOT is required for BCalTagEfficiency
LINK_LIBRARIES( ${ROOT_MATHMORE_LIBRARY} )

FIND_PACKAGE( FORTRAN_LIBRARY REQUIRED )
LINK_LIBRARIES( ${FORTRAN_LIBRARY} )
GET_FILENAME_COMPONENT( FORTRAN_LIBRARY_FULL_PATH ${FORTRAN_LIBRARY} REALPATH ) 
INSTALL( PROGRAMS ${FORTRAN_LIBRARY_FULL_PATH} DESTINATION lib ) # PROGRAMS is like FILES but sets executable permissions


# -------- compiler flags ---------------------------------------------------
IF( FORTRAN_COMPILER_NAME STREQUAL "g77" )
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-globals" ) # -ff2c ?
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dg77Fortran" ) # -Df2cFortran
ELSE() # gfortran / gfortran-32
    SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DgFortran" )
ENDIF()
SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DNDEBUG" )
SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -DNDEBUG" )
# ---------------------------------------------------------------------------



### DOCUMENTATION ############################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )
IF( INSTALL_DOC )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY ##################################################################

ADD_SUBDIRECTORY( ./src )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )

