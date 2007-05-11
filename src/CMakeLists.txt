######################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
######################################

####################################################################################################
### INCLUDES
####################################################################################################

INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Analysis/EventShapes/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Analysis/SatoruJetFinder/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Calibration/AbsCalibration/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/CaloDigi/LDCCaloDigi/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Clustering/ClusterCheater/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Clustering/ClusterCheater5_3/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Clustering/NNClustering/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Clustering/TrackwiseClustering/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Clustering/PhotonFinderKit/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Pflow/Wolf/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Pflow/TrackBasedPFlow/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/TrackDigi/FTDDigi/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/TrackDigi/TPCDigi/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/TrackDigi/VTXDigi/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Tracking/BrahmsTracking/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Tracking/FullLDCTracking/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Tracking/SiliconTracking/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Tracking/TrackCheater/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/Tracking/VertexTracking/include" )


####################################################################################################
### SOURCES
####################################################################################################

####################################################################################################
# EventShapes

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/EventShapes/src EventShapes_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${EventShapes_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# SatoruJetFinder

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/SatoruJetFinder/src SatoruJetFinder_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${SatoruJetFinder_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )

# fortran sources
FILE( GLOB SatoruJetFinder_f77_srcs ../Analysis/SatoruJetFinder/src/fortsrc/*.f )

# fortran flags
SET_SOURCE_FILES_PROPERTIES( ${f77_srcs} PROPERTIES
	COMPILE_FLAGS "-Wno-globals -ff2c" )


####################################################################################################
# AbsCalibration

# cxx sources
AUX_SOURCE_DIRECTORY( ../Calibration/AbsCalibration/src AbsCalibration_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${AbsCalibration_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# LDCCaloDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../CaloDigi/LDCCaloDigi/src LDCCaloDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${LDCCaloDigi_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# ClusterCheater

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/ClusterCheater/src  ClusterCheater_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ClusterCheater_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# ClusterCheater5_3

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/ClusterCheater5_3/src  ClusterCheater5_3_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ClusterCheater5_3_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# NNClustering

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/NNClustering/src NNClustering_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${NNClustering_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# TrackwiseClustering

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/TrackwiseClustering/src TrackwiseClustering_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackwiseClustering_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# PhotonFinderKit

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/PhotonFinderKit/src PhotonFinderKit_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${PhotonFinderKit_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# Wolf

# cxx sources
AUX_SOURCE_DIRECTORY( ../Pflow/Wolf/src Wolf_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${Wolf_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# TrackBasedPFlow

# cxx sources
AUX_SOURCE_DIRECTORY( ../Pflow/TrackBasedPFlow/src TrackBasedPFlow_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackBasedPFlow_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# FTDDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/FTDDigi/src FTDDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${FTDDigi_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )



####################################################################################################
# TPCDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/TPCDigi/src TPCDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TPCDigi_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# VTXDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/VTXDigi/src VTXDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${VTXDigi_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# BrahmsTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/BrahmsTracking/src BrahmsTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${BrahmsTracking_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic"
	# required for the fortran interface to work
	COMPILE_FLAGS "-Df2cFortran" )

# c sources
AUX_SOURCE_DIRECTORY( ../Tracking/BrahmsTracking/src/C BrahmsTracking_c_srcs )

SET_SOURCE_FILES_PROPERTIES( ${BrahmsTracking_c_srcs} PROPERTIES
	# require proper c
	COMPILE_FLAGS "-Wall -ansi -pedantic"
	COMPILE_FLAGS "-Wno-globals" )

# fortran sources
FILE( GLOB BrahmsTracking_f77_srcs ../Tracking/BrahmsTracking/src/f77/*.F )

SET_SOURCE_FILES_PROPERTIES( ${BrahmsTracking_f77_srcs} PROPERTIES
	COMPILE_FLAGS "-Wno-globals -ff2c" )


####################################################################################################
# FullLDCTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/FullLDCTracking/src FullLDCTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${FullLDCTracking_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# SiliconTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/SiliconTracking/src SiliconTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${SiliconTracking_cxx_srcs} PROPERTIES
	# FIXME require proper c++
	#COMPILE_FLAGS "-Wall -ansi -pedantic" )
	COMPILE_FLAGS "-Wall -ansi" )


####################################################################################################
# TrackCheater

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/TrackCheater/src TrackCheater_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackCheater_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


####################################################################################################
# VertexTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/VertexTracking/src VertexTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${VertexTracking_cxx_srcs} PROPERTIES
	# require proper c++
	COMPILE_FLAGS "-Wall -ansi -pedantic" )


# set library sources
SET( lib_srcs
	${EventShapes_cxx_srcs}
	${SatoruJetFinder_cxx_srcs}
	${SatoruJetFinder_f77_srcs}
	${AbsCalibration_cxx_srcs}
	${LDCCaloDigi_cxx_srcs}
	${ClusterCheater_cxx_srcs}
	${ClusterCheater5_3_cxx_srcs}
	${NNClustering_cxx_srcs}
	${TrackwiseClustering_cxx_srcs}
	${PhotonFinderKit_cxx_srcs}
	${Wolf_cxx_srcs}
	${TrackBasedPFlow_cxx_srcs}
	${FTDDigi_cxx_srcs}
	${TPCDigi_cxx_srcs}
	${VTXDigi_cxx_srcs}
	${BrahmsTracking_cxx_srcs}
	${BrahmsTracking_c_srcs}
	${BrahmsTracking_f77_srcs}
	${FullLDCTracking_cxx_srcs}
	${SiliconTracking_cxx_srcs}
	${TrackCheater_cxx_srcs}
	${VertexTracking_cxx_srcs}
)

####################################################################################################
### LIBRARY
####################################################################################################

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# workaround for fortran cmake bug
ADD_CUSTOM_TARGET( fix_f77_cmake_bug

	# SatoruJetFinder fix
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > original.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > this.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > YCUT.mod.proxy

	# BrahmsTracking fix
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > Monte.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > data.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > a.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > A.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > GEAR.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > of.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > helix.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > as.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > GRNDM.mod.proxy
	COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > ranecu.mod.proxy

	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	COMMENT "Fixing cmake fortran bug..." VERBATIM )

# establish dependency for fixing fortran bug
ADD_DEPENDENCIES( lib_${PROJECT_NAME} fix_f77_cmake_bug )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
	# create *nix style library versions + symbolic links
	VERSION ${${PROJECT_NAME}_VERSION}
	SOVERSION ${${PROJECT_NAME}_SOVERSION}
	# allow creating static and shared libs without conflicts
	CLEAN_DIRECT_OUTPUT 1
	# linker must be explicitly set due to mixed sources
	LINKER_LANGUAGE CXX
	# avoid conflicts between library and binary target names
	OUTPUT_NAME ${PROJECT_NAME}	)

# FIXME ******************************************************************
# do not store find results in cache
SET( G2C_LIB G2C_LIB-NOTFOUND )
FIND_LIBRARY( G2C_LIB
    NAMES g2c
    PATHS /usr/lib/gcc-lib/i386-redhat-linux/3.2.3
)
IF( NOT G2C_LIB )
    MESSAGE( FATAL_ERROR "failed to find g2c library!!" )
ENDIF()
MARK_AS_ADVANCED( G2C_LIB )

# link MarlinReco library to lib g2c
TARGET_LINK_LIBRARIES( lib_${PROJECT_NAME} ${G2C_LIB} )
# FIXME END ***************************************************************

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
		OWNER_READ OWNER_WRITE OWNER_EXECUTE
		GROUP_READ GROUP_EXECUTE
		WORLD_READ WORLD_EXECUTE )

####################################################################################################
### DEPENDENCIES
####################################################################################################

# load macro
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
	MESSAGE( STATUS "-------------------------------------------------------------------------------" )
	MESSAGE( STATUS "Change a module with: cmake -D<ModuleName>_HOME=<Path_to_Module>" )
	MESSAGE( STATUS )

	# load macro
	IF( NOT EXISTS "${CMAKE_MODULE_PATH}/LoadPackageMacro.cmake" )
		MESSAGE( FATAL_ERROR "\nSorry, could not find LoadPackageMacro.cmake...\n" 
			"Please set CMAKE_MODULE_PATH correctly with: "
			"cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
	ENDIF()
	INCLUDE( "${CMAKE_MODULE_PATH}/LoadPackageMacro.cmake" )
ENDIF()

# project dependencies
IF( DEFINED ${PROJECT_NAME}_DEPENDS )
	SEPARATE_ARGUMENTS( ${PROJECT_NAME}_DEPENDS )
	MARK_AS_ADVANCED( ${PROJECT_NAME}_DEPENDS )
	FOREACH( req_pkg ${${PROJECT_NAME}_DEPENDS} )
		LOAD_PACKAGE( ${req_pkg} REQUIRED )
	ENDFOREACH()
ENDIF()

# user defined dependencies
IF( DEFINED BUILD_WITH )
	SEPARATE_ARGUMENTS( BUILD_WITH )
	MARK_AS_ADVANCED( BUILD_WITH )
	FOREACH( opt_pkg ${BUILD_WITH} )
		LOAD_PACKAGE( ${opt_pkg} )
	ENDFOREACH()
	SET( BUILD_WITH "${BUILD_WITH}" CACHE STRING
		"Build ${PROJECT_NAME} with these optional packages" FORCE )
ENDIF()

# user defined dependencies
IF( DEFINED LINK_WITH )
	SEPARATE_ARGUMENTS( LINK_WITH )
	MARK_AS_ADVANCED( LINK_WITH )
	FOREACH( lnk_pkg ${LINK_WITH} )
		LOAD_PACKAGE( ${lnk_pkg} LINK_ONLY )
	ENDFOREACH()
	SET( LINK_WITH "${LINK_WITH}" CACHE STRING
		"Link ${PROJECT_NAME} with these optional packages" FORCE )
ENDIF()

IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
	MESSAGE( STATUS "-------------------------------------------------------------------------------" )
ENDIF()

