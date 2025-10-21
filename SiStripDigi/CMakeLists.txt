#############################################################################
#
# CMAKE file for building Marlin package - SiStripDigi
#
# @author Zbynek Drasal
# @author J. Duarte Campderros (updates)
#
#############################################################################
# CMake compatibility issues: don't modify this, please!
CMAKE_MINIMUM_REQUIRED( VERSION 2.6.2 FATAL_ERROR )

# Project name
PROJECT( SiStripDigi )


# Project version
SET( ${PROJECT_NAME}_MAJOR_VERSION 2 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )


SET( ILCUTIL_DIR /usr/local/ilcsoft/v01-08/ilcutil/trunk )
SET( GEAR_DIR /usr/local/ilcsoft/v01-08/gear/trunk )
SET( LCIO_DIR /usr/local/ilcsoft/v01-08/lcio/trunk )


#####
ADD_DEFINITIONS( -DROOT_OUTPUT_LAND -DROOT_OUTPUT)
ADD_DEFINITIONS( -DUSE_GEAR )


### DEPENDENCIES #######################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES ILCTEST REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

# required packages
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( GEAR REQUIRED )
FIND_PACKAGE( Marlin REQUIRED )

# export SiStripDigi_DEPENDS_INCLUDE_DIRS to SiStripDigiConfig.cmake
SET( SiStripDigi_DEPENDS_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS} ${Marlin_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS} ${GEAR_INCLUDE_DIRS} ${streamlog_INCLUDE_DIRS} )
SET( SiStripDigi_DEPENDS_LIBRARY_DIRS ${ROOT_LIBRARY_DIRS} ${CLHEP_LIBRARY_DIRS} ${Marlin_LIBRARY_DIRS} ${LCIO_LIBRARY_DIRS} ${GEAR_LIBRARY_DIRS} ${streamlog_LIBRARY_DIRS} )
SET( SiStripDigi_DEPENDS_LIBRARIES ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} ${Marlin_LIBRARIES} ${LCIO_LIBRARIES} ${GEAR_LIBRARIES} ${streamlog_LIBRARIES} )

INCLUDE_DIRECTORIES( ${SiStripDigi_DEPENDS_INCLUDE_DIRS} )
LINK_LIBRARIES( ${SiStripDigi_DEPENDS_LIBRARIES} )


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )

    FIND_PACKAGE( Doxygen )

    IF( DOXYGEN_FOUND )
        ADD_SUBDIRECTORY( ./doc )
    ELSE()
        MESSAGE( SEND_ERROR "Could not find doxygen required to build documentation" )
        MESSAGE( "Please install doxygen or set INSTALL_DOC to OFF" )
        MESSAGE( "" )
    ENDIF()

ENDIF()


### LIBRARY AND TOOLS #######################################################

ADD_SUBDIRECTORY( ./source )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( SiStripDigiConfig.cmake SiStripDigiConfigVersion.cmake SiStripDigiLibDeps.cmake )



