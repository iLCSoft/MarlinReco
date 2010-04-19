###############################################
# cmake configuration file for PrepareRECParticles
# @author Jan Engels, DESY
###############################################

SET( PrepareRECParticles_FOUND FALSE )
MARK_AS_ADVANCED( PrepareRECParticles_FOUND )

# do not store find results in cache
SET( PrepareRECParticles_INCLUDE_DIR PrepareRECParticles_INCLUDE_DIR-NOTFOUND )

FIND_PATH( PrepareRECParticles_INCLUDE_DIR
	NAMES PrepareRECParticles.h
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/PrepareRECParticles
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT PrepareRECParticles_INCLUDE_DIR )
    MESSAGE( STATUS "Check for PrepareRECParticles: ${PrepareRECParticles_HOME}"
					" -- failed to find PrepareRECParticles include directory!!" )
ELSE( NOT PrepareRECParticles_INCLUDE_DIR )
    MARK_AS_ADVANCED( PrepareRECParticles_INCLUDE_DIR )
ENDIF( NOT PrepareRECParticles_INCLUDE_DIR )


# do not store find results in cache
SET( PrepareRECParticles_LIB PrepareRECParticles_LIB-NOTFOUND )

FIND_LIBRARY( PrepareRECParticles_LIB
	NAMES PrepareRECParticles
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/PrepareRECParticles
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT PrepareRECParticles_LIB )
    MESSAGE( STATUS "Check for PrepareRECParticles: ${PrepareRECParticles_HOME}"
					" -- failed to find PrepareRECParticles library!!" )
ELSE( NOT PrepareRECParticles_LIB )
    MARK_AS_ADVANCED( PrepareRECParticles_LIB )
ENDIF( NOT PrepareRECParticles_LIB )


# set variables and display results
IF( PrepareRECParticles_INCLUDE_DIR AND PrepareRECParticles_LIB )
    SET( PrepareRECParticles_FOUND TRUE )
    SET( PrepareRECParticles_INCLUDE_DIRS ${PrepareRECParticles_INCLUDE_DIR} )
	SET( PrepareRECParticles_LIBRARIES ${PrepareRECParticles_LIB} )
	MARK_AS_ADVANCED( PrepareRECParticles_INCLUDE_DIRS PrepareRECParticles_LIBRARIES )
	MESSAGE( STATUS "Check for PrepareRECParticles: ${PrepareRECParticles_HOME} -- works" )
ELSE( PrepareRECParticles_INCLUDE_DIR AND PrepareRECParticles_LIB )
	IF( PrepareRECParticles_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for PrepareRECParticles: ${PrepareRECParticles_HOME} -- failed!!" )
    ELSE( PrepareRECParticles_FIND_REQUIRED )
        MESSAGE( STATUS "Check for PrepareRECParticles: ${PrepareRECParticles_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( PrepareRECParticles_FIND_REQUIRED )
ENDIF( PrepareRECParticles_INCLUDE_DIR AND PrepareRECParticles_LIB )
