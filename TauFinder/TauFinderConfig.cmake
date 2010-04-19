###############################################
# cmake configuration file for TauFinder
# @author Jan Engels, DESY
###############################################

SET( TauFinder_FOUND FALSE )
MARK_AS_ADVANCED( TauFinder_FOUND )

# do not store find results in cache
SET( TauFinder_INCLUDE_DIR TauFinder_INCLUDE_DIR-NOTFOUND )

FIND_PATH( TauFinder_INCLUDE_DIR
	NAMES TauFinder.h
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/TauFinder
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT TauFinder_INCLUDE_DIR )
    MESSAGE( STATUS "Check for TauFinder: ${TauFinder_HOME}"
					" -- failed to find TauFinder include directory!!" )
ELSE( NOT TauFinder_INCLUDE_DIR )
    MARK_AS_ADVANCED( TauFinder_INCLUDE_DIR )
ENDIF( NOT TauFinder_INCLUDE_DIR )


# do not store find results in cache
SET( TauFinder_LIB TauFinder_LIB-NOTFOUND )

FIND_LIBRARY( TauFinder_LIB
	NAMES TauFinder
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/TauFinder
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT TauFinder_LIB )
    MESSAGE( STATUS "Check for TauFinder: ${TauFinder_HOME}"
					" -- failed to find TauFinder library!!" )
ELSE( NOT TauFinder_LIB )
    MARK_AS_ADVANCED( TauFinder_LIB )
ENDIF( NOT TauFinder_LIB )


# set variables and display results
IF( TauFinder_INCLUDE_DIR AND TauFinder_LIB )
    SET( TauFinder_FOUND TRUE )
    SET( TauFinder_INCLUDE_DIRS ${TauFinder_INCLUDE_DIR} )
	SET( TauFinder_LIBRARIES ${TauFinder_LIB} )
	MARK_AS_ADVANCED( TauFinder_INCLUDE_DIRS TauFinder_LIBRARIES )
	MESSAGE( STATUS "Check for TauFinder: ${TauFinder_HOME} -- works" )
ELSE( TauFinder_INCLUDE_DIR AND TauFinder_LIB )
	IF( TauFinder_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for TauFinder: ${TauFinder_HOME} -- failed!!" )
    ELSE( TauFinder_FIND_REQUIRED )
        MESSAGE( STATUS "Check for TauFinder: ${TauFinder_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( TauFinder_FIND_REQUIRED )
ENDIF( TauFinder_INCLUDE_DIR AND TauFinder_LIB )
