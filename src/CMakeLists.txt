######################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
######################################


# MarlinReco Packages
SET( MarlinReco_PKGS
    "Analysis/CLICPfoSelector"
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
    "Tracking/KinkFinder" 
    "Tracking/V0Finder"
)

IF( CERNLIB_FOUND ) 
SET( MarlinReco_PKGS ${MarlinReco_PKGS}
    "Tracking/BrahmsTracking"
    "Tracking/FullLDCTracking"
    "Tracking/SiliconTracking"
    "Analysis/SatoruJetFinder"
    "Tracking/TrackCheater"
    "Analysis/BCalTagEfficiency"
    "Analysis/MarlinKinfit"
    "PFOID"
    "Analysis/EventShapes"
)
ELSE()
  MESSAGE( WARNING "=====================================================================================
  MarlinReco will be built without cernlib - the following packages will not be available : 
    Tracking/BrahmsTracking
    Tracking/FullLDCTracking
    Tracking/SiliconTracking
    Analysis/SatoruJetFinder
    Tracking/TrackCheater
    Analysis/BCalTagEfficiency
    Analysis/MarlinKinfit
    PFOID
    Analysis/EventShapes
 In order to build these packages specify MARLIN_CERNLIB=on and CERNLIB_DIR !! 
 =============================================================================================  " 
)

ENDIF(CERNLIB_FOUND)

####################################################################################################
### INCLUDES and SOURCES
####################################################################################################

FOREACH( pkg ${MarlinReco_PKGS} )
    MESSAGE( STATUS "Including MarlinReco package: ${pkg}" )
    INCLUDE_DIRECTORIES( BEFORE "${PROJECT_SOURCE_DIR}/${pkg}/include" )
    INSTALL_DIRECTORY( "${PROJECT_SOURCE_DIR}/${pkg}/include" DESTINATION . )

    AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/${pkg}/src MarlinReco_cxx_srcs )

ENDFOREACH()

IF( CERNLIB_FOUND ) 
# C sources
AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/Tracking/BrahmsTracking/src/C MarlinReco_c_srcs )

# Fortran sources
 FILE( GLOB SatoruJetFinder_f77_srcs ${PROJECT_SOURCE_DIR}/Analysis/SatoruJetFinder/src/fortsrc/*.f )
 FILE( GLOB BrahmsTracking_f77_srcs ${PROJECT_SOURCE_DIR}/Tracking/BrahmsTracking/src/f77/*.F )
 LIST( APPEND MarlinReco_f77_srcs "${PROJECT_SOURCE_DIR}/Analysis/BCalTagEfficiency/Fortran/track-to-bcal.F" )
 LIST( APPEND MarlinReco_f77_srcs ${SatoruJetFinder_f77_srcs} )
 LIST( APPEND MarlinReco_f77_srcs ${BrahmsTracking_f77_srcs} )
ENDIF(CERNLIB_FOUND)


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
ADD_SHARED_LIBRARY( MarlinReco ${MarlinReco_cxx_srcs} ${MarlinReco_c_srcs} ${MarlinReco_f77_srcs} )
INSTALL_SHARED_LIBRARY( MarlinReco DESTINATION lib )

