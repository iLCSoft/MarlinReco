########################################################
# cmake file for building MarlinReco
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED( VERSION 2.6.2 FATAL_ERROR )
########################################################

# fix for lcio: all c++ sources need to be compiled with -Wno-long-long
SET( ENV{CXXFLAGS} "-Wno-long-long $ENV{CXXFLAGS}" )

# require proper c and c++
# -pedantic cannot be set globally because it cannot be removed for certain
# source files only: grep pedantic src/CMakeLists.txt
SET( ENV{CXXFLAGS} "-Wall -ansi $ENV{CXXFLAGS}" ) # FIXME -pedantic
SET( ENV{CFLAGS} "-w -ansi -pedantic $ENV{CFLAGS}" ) # -w suppress all warning messages


####################################################################################################

OPTION( MARLINRECO_FORTRAN "Set to ON to build MarlinReco fortran sources" ON )

# project name
IF( MARLINRECO_FORTRAN )
    PROJECT( MarlinReco C CXX Fortran )
ELSE()
    PROJECT( MarlinReco C CXX )
ENDIF()


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 3 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.4 REQUIRED )
FIND_PACKAGE( GSL REQUIRED )
FIND_PACKAGE( ROOT 5.27 REQUIRED COMPONENTS MathMore )

OPTION( MARLINRECO_AIDA "Set to ON to build MarlinReco with AIDA" ON )
IF( MARLINRECO_AIDA )
    FIND_PACKAGE( AIDA )
    IF( NOT AIDA_FOUND )
        MESSAGE( SEND_ERROR "AIDA not found. Please set MARLINRECO_AIDA to OFF or AIDA_DIR=/path/to/AIDA" )
    ENDIF()
ENDIF()

FOREACH( pkg Marlin MarlinUtil GSL AIDA ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

# MathMore library from ROOT is required for BCalTagEfficiency
LINK_LIBRARIES( ${ROOT_MATHMORE_LIBRARY} )

IF( MARLINRECO_FORTRAN )

    FIND_PACKAGE( CERNLIB REQUIRED )
    INCLUDE_DIRECTORIES( ${CERNLIB_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${CERNLIB_LIBRARIES} )

    FIND_PACKAGE( FORTRAN_LIBRARY REQUIRED )
    LINK_LIBRARIES( ${FORTRAN_LIBRARY} )
    GET_FILENAME_COMPONENT( FORTRAN_LIBRARY_FULL_PATH ${FORTRAN_LIBRARY} REALPATH ) 
    INSTALL( PROGRAMS ${FORTRAN_LIBRARY_FULL_PATH} DESTINATION lib ) # PROGRAMS is like FILES but sets executable permissions


    # -------- compiler flags ---------------------------------------------------
    IF( FORTRAN_COMPILER_NAME STREQUAL "g77" )
        SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-globals" ) # -ff2c ?
        #SET( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffortran-bounds-check" )
        SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dg77Fortran" ) # -Df2cFortran
    ELSE() # gfortran / gfortran-32
        SET( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )
        #SET( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check" )
        SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DgFortran" )
    ENDIF()
    SET( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DNDEBUG" )
    SET( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL} -DNDEBUG" )
    # ---------------------------------------------------------------------------

ENDIF( MARLINRECO_FORTRAN )



### DOCUMENTATION ############################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )
IF( INSTALL_DOC )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY ##################################################################

MACRO( ADD_MARLINRECO_PKG _input_dir )
    MESSAGE( STATUS "Including MarlinReco package: ${_input_dir}" )
    INCLUDE_DIRECTORIES( BEFORE ${_input_dir}/include )
    AUX_SOURCE_DIRECTORY( ${_input_dir}/src MarlinReco_cxx_srcs )
    #INSTALL_DIRECTORY( "${_input_dir}/include" DESTINATION . FILES_MATCHING PATTERN "*.h" )
    #INSTALL_DIRECTORY( ${_input_dir}/include/ DESTINATION include/marlinreco FILES_MATCHING PATTERN "*.h" )
ENDMACRO()


ADD_MARLINRECO_PKG( ./Analysis/IsolatedLeptonFinder )
ADD_MARLINRECO_PKG( ./Analysis/CLICPfoSelector )
ADD_MARLINRECO_PKG( ./Analysis/RecoMCTruthLink )
ADD_MARLINRECO_PKG( ./Analysis/ZFinder )
ADD_MARLINRECO_PKG( ./Calibration/AbsCalibration )
ADD_MARLINRECO_PKG( ./CaloDigi/LDCCaloDigi )
ADD_MARLINRECO_PKG( ./Clustering/BCalReco) 
ADD_MARLINRECO_PKG( ./Clustering/ClusterCheater )
ADD_MARLINRECO_PKG( ./Clustering/ClusterCheater5_3 )
ADD_MARLINRECO_PKG( ./Clustering/NNClustering )
#ADD_MARLINRECO_PKG( ./Clustering/TrackwiseClustering )
ADD_MARLINRECO_PKG( ./Clustering/PhotonFinderKit )
ADD_MARLINRECO_PKG( ./Clustering/hybridEcalSplitter )
#ADD_MARLINRECO_PKG( ./Pflow/Wolf )
#ADD_MARLINRECO_PKG( ./Pflow/TrackBasedPFlow )
ADD_MARLINRECO_PKG( ./TrackDigi/TPCDigi )
ADD_MARLINRECO_PKG( ./TrackDigi/FTDDigi )
ADD_MARLINRECO_PKG( ./TrackDigi/ETDDigi )
ADD_MARLINRECO_PKG( ./TrackDigi/VTXDigi )
ADD_MARLINRECO_PKG( ./TrackDigi/FPCCDDigi )
ADD_MARLINRECO_PKG( ./Tracking/KinkFinder )
ADD_MARLINRECO_PKG( ./Tracking/V0Finder )
ADD_MARLINRECO_PKG( ./PFOID )


IF( MARLINRECO_FORTRAN )
    ADD_MARLINRECO_PKG( ./Analysis/BCalTagEfficiency )
    ADD_MARLINRECO_PKG( ./Analysis/EventShapes )
    ADD_MARLINRECO_PKG( ./Analysis/SatoruJetFinder )
    ADD_MARLINRECO_PKG( ./Tracking/BrahmsTracking )
    ADD_MARLINRECO_PKG( ./Tracking/FullLDCTracking )
    ADD_MARLINRECO_PKG( ./Tracking/SiliconTracking )
    ADD_MARLINRECO_PKG( ./Tracking/TrackCheater )
    AUX_SOURCE_DIRECTORY( ./Tracking/BrahmsTracking/src/C MarlinReco_c_srcs )
    FILE( GLOB SatoruJetFinder_f77_srcs ./Analysis/SatoruJetFinder/src/fortsrc/*.f )
    FILE( GLOB BrahmsTracking_f77_srcs ./Tracking/BrahmsTracking/src/f77/*.F )
    LIST( APPEND MarlinReco_f77_srcs ./Analysis/BCalTagEfficiency/Fortran/track-to-bcal.F )
    LIST( APPEND MarlinReco_f77_srcs ${SatoruJetFinder_f77_srcs} )
    LIST( APPEND MarlinReco_f77_srcs ${BrahmsTracking_f77_srcs} )
ENDIF()


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


# create library
ADD_SHARED_LIBRARY( MarlinReco ${MarlinReco_cxx_srcs} ${MarlinReco_c_srcs} ${MarlinReco_f77_srcs} )
INSTALL_SHARED_LIBRARY( MarlinReco DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

