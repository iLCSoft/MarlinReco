######################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
######################################


# MarlinReco Packages
SET( MarlinReco_PKGS
    "Analysis/EventShapes"
    "Analysis/BCalTagEfficiency"
    "Analysis/SatoruJetFinder"
    "Analysis/MarlinKinfit"
    "Analysis/RecoMCTruthLink"
    "Analysis/ZFinder"
    "Calibration/AbsCalibration"
    "CaloDigi/LDCCaloDigi"
    "Clustering/ClusterCheater" 
    "Clustering/ClusterCheater5_3"
    "Clustering/NNClustering"
    "Clustering/TrackwiseClustering"
    "Clustering/PhotonFinderKit"
    "Pflow/Wolf"
    "Pflow/TrackBasedPFlow"
    "TrackDigi/TPCDigi"
    "TrackDigi/FTDDigi"
    "TrackDigi/ETDDigi"
    "TrackDigi/VTXDigi"
    "TrackDigi/FPCCDDigi"
    "Tracking/BrahmsTracking"
    "Tracking/FullLDCTracking"
    "Tracking/SiliconTracking"
    "Tracking/TrackCheater"
    "Tracking/KinkFinder" 
    "Tracking/V0Finder"
    "PFOID"
)



####################################################################################################
### INCLUDES and SOURCES
####################################################################################################

FOREACH( pkg ${MarlinReco_PKGS} )
    MESSAGE( STATUS "Including MarlinReco package: ${pkg}" )
    INCLUDE_DIRECTORIES( BEFORE "${PROJECT_SOURCE_DIR}/${pkg}/include" )
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/${pkg}/include" DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE
        PATTERN "*/.svn*" EXCLUDE
    )

    AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/${pkg}/src MarlinReco_cxx_srcs )

ENDFOREACH()

# C sources
AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/Tracking/BrahmsTracking/src/C MarlinReco_c_srcs )

# Fortran sources
FILE( GLOB SatoruJetFinder_f77_srcs ${PROJECT_SOURCE_DIR}/Analysis/SatoruJetFinder/src/fortsrc/*.f )
FILE( GLOB BrahmsTracking_f77_srcs ${PROJECT_SOURCE_DIR}/Tracking/BrahmsTracking/src/f77/*.F )
LIST( APPEND MarlinReco_f77_srcs "${PROJECT_SOURCE_DIR}/Analysis/BCalTagEfficiency/Fortran/track-to-bcal.F" )
LIST( APPEND MarlinReco_f77_srcs ${SatoruJetFinder_f77_srcs} )
LIST( APPEND MarlinReco_f77_srcs ${BrahmsTracking_f77_srcs} )


#STRING( REPLACE "-pedantic" "" CXX_FLAGS_NOT_PEDANTIC ${CMAKE_CXX_FLAGS} )
#MESSAGE( STATUS "CXX_FLAGS_NOT_PEDANTIC: ${CXX_FLAGS_NOT_PEDANTIC}" )
SET_SOURCE_FILES_PROPERTIES( ${MarlinReco_cxx_srcs} PROPERTIES COMPILE_FLAGS "-pedantic" )

# FIXME sources which do not compile with -pedantic
SET_SOURCE_FILES_PROPERTIES(
    "${PROJECT_SOURCE_DIR}/Clustering/PhotonFinderKit/src/EMShowerFinder.cc"
    "${PROJECT_SOURCE_DIR}/Clustering/PhotonFinderKit/src/KITutil.cc"
    "${PROJECT_SOURCE_DIR}/Clustering/PhotonFinderKit/src/KIT.cc"
    "${PROJECT_SOURCE_DIR}/Pflow/TrackBasedPFlow/src/TrackBasedPFlow.cc"
    "${PROJECT_SOURCE_DIR}/Clustering/NNClustering/src/NNClusterProcessor.cc"
    "${PROJECT_SOURCE_DIR}/Tracking/SiliconTracking/src/SiliconTracking.cc"
    "${PROJECT_SOURCE_DIR}/Tracking/KinkFinder/src/KinkFinder.cc"
    PROPERTIES COMPILE_FLAGS ""
)

####################################################################################################
### LIBRARY
####################################################################################################

#MESSAGE( STATUS "MarlinReco_cxx_srcs: ${MarlinReco_cxx_srcs}" )

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${MarlinReco_cxx_srcs} ${MarlinReco_c_srcs} ${MarlinReco_f77_srcs} )

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
    OUTPUT_NAME ${PROJECT_NAME} )

# link MarlinReco library to libg2c/libgfortran
TARGET_LINK_LIBRARIES( lib_${PROJECT_NAME} ${FORTRAN_LIB} )
IF(APPLE)
  TARGET_LINK_LIBRARIES( lib_${PROJECT_NAME} streamlog )
ENDIF(APPLE)

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

