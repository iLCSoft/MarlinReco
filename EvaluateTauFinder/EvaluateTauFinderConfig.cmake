###############################################
# cmake configuration file for EvaluateTauFinder
# @author Jan Engels, DESY
###############################################

SET( EvaluateTauFinder_FOUND FALSE )
MARK_AS_ADVANCED( EvaluateTauFinder_FOUND )

# do not store find results in cache
SET( EvaluateTauFinder_INCLUDE_DIR EvaluateTauFinder_INCLUDE_DIR-NOTFOUND )

FIND_PATH( EvaluateTauFinder_INCLUDE_DIR
	NAMES EvaluateTauFinder.h
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/EvaluateTauFinder
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT EvaluateTauFinder_INCLUDE_DIR )
    MESSAGE( STATUS "Check for EvaluateTauFinder: ${EvaluateTauFinder_HOME}"
					" -- failed to find EvaluateTauFinder include directory!!" )
ELSE( NOT EvaluateTauFinder_INCLUDE_DIR )
    MARK_AS_ADVANCED( EvaluateTauFinder_INCLUDE_DIR )
ENDIF( NOT EvaluateTauFinder_INCLUDE_DIR )


# do not store find results in cache
SET( EvaluateTauFinder_LIB EvaluateTauFinder_LIB-NOTFOUND )

FIND_LIBRARY( EvaluateTauFinder_LIB
	NAMES EvaluateTauFinder
	PATHS /afs/cern.ch/eng/clic/work/amunnich/myMarlin/processors/EvaluateTauFinder
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT EvaluateTauFinder_LIB )
    MESSAGE( STATUS "Check for EvaluateTauFinder: ${EvaluateTauFinder_HOME}"
					" -- failed to find EvaluateTauFinder library!!" )
ELSE( NOT EvaluateTauFinder_LIB )
    MARK_AS_ADVANCED( EvaluateTauFinder_LIB )
ENDIF( NOT EvaluateTauFinder_LIB )


# set variables and display results
IF( EvaluateTauFinder_INCLUDE_DIR AND EvaluateTauFinder_LIB )
    SET( EvaluateTauFinder_FOUND TRUE )
    SET( EvaluateTauFinder_INCLUDE_DIRS ${EvaluateTauFinder_INCLUDE_DIR} )
	SET( EvaluateTauFinder_LIBRARIES ${EvaluateTauFinder_LIB} )
	MARK_AS_ADVANCED( EvaluateTauFinder_INCLUDE_DIRS EvaluateTauFinder_LIBRARIES )
	MESSAGE( STATUS "Check for EvaluateTauFinder: ${EvaluateTauFinder_HOME} -- works" )
ELSE( EvaluateTauFinder_INCLUDE_DIR AND EvaluateTauFinder_LIB )
	IF( EvaluateTauFinder_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for EvaluateTauFinder: ${EvaluateTauFinder_HOME} -- failed!!" )
    ELSE( EvaluateTauFinder_FIND_REQUIRED )
        MESSAGE( STATUS "Check for EvaluateTauFinder: ${EvaluateTauFinder_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( EvaluateTauFinder_FIND_REQUIRED )
ENDIF( EvaluateTauFinder_INCLUDE_DIR AND EvaluateTauFinder_LIB )
