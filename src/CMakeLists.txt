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
    "Tracking/BrahmsTracking"
    "Tracking/FullLDCTracking"
    "Tracking/SiliconTracking"
    "Tracking/TrackCheater"
    "Tracking/KinkFinder" 
    "Tracking/V0Finder"
    "PFOID"
)



####################################################################################################
### INCLUDES
####################################################################################################

FOREACH( pkg ${MarlinReco_PKGS} )
    MESSAGE( STATUS "Including MarlinReco package: ${pkg}" )
    INCLUDE_DIRECTORIES( "../${pkg}/include" )
    INSTALL( DIRECTORY "../${pkg}/include" DESTINATION .
        PATTERN "*~" EXCLUDE
        PATTERN "*CVS*" EXCLUDE
        PATTERN "*/.svn*" EXCLUDE
        )
ENDFOREACH()



####################################################################################################
### SOURCES
####################################################################################################



####################################################################################################
# BCalTagEfficiency

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/BCalTagEfficiency/src BCalTagEfficiency_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${BCalTagEfficiency_cxx_srcs} PROPERTIES
    COMPILE_FLAGS "-Wall" )

SET( BCalTagEfficiency_f77_srcs ../Analysis/BCalTagEfficiency/Fortran/track-to-bcal.F )

SET( BCalTagEfficiency_srcs ${BCalTagEfficiency_cxx_srcs} ${BCalTagEfficiency_f77_srcs} )

####################################################################################################
# EventShapes

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/EventShapes/src EventShapes_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${EventShapes_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( EventShapes_srcs ${EventShapes_cxx_srcs} )

####################################################################################################
# MarlinKinFit

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/MarlinKinfit/src MarlinKinfit_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${MarlinKinfit_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( MarlinKinfit_srcs ${MarlinKinfit_cxx_srcs} )


####################################################################################################
# RecoMCTruthLink

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/RecoMCTruthLink/src RecoMCTruthLink_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${RecoMCTruthLink_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( RecoMCTruthLink_srcs ${RecoMCTruthLink_cxx_srcs} )


####################################################################################################
# SatoruJetFinder

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/SatoruJetFinder/src SatoruJetFinder_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${SatoruJetFinder_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

# fortran sources
FILE( GLOB SatoruJetFinder_f77_srcs ../Analysis/SatoruJetFinder/src/fortsrc/*.f )

SET( SatoruJetFinder_srcs ${SatoruJetFinder_cxx_srcs} ${SatoruJetFinder_f77_srcs} )


####################################################################################################
# ZFinder

# cxx sources
AUX_SOURCE_DIRECTORY( ../Analysis/ZFinder/src ZFinder_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ZFinder_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( ZFinder_srcs ${ZFinder_cxx_srcs} )


####################################################################################################
# AbsCalibration

# cxx sources
AUX_SOURCE_DIRECTORY( ../Calibration/AbsCalibration/src AbsCalibration_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${AbsCalibration_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( AbsCalibration_srcs ${AbsCalibration_cxx_srcs} )


####################################################################################################
# LDCCaloDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../CaloDigi/LDCCaloDigi/src LDCCaloDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${LDCCaloDigi_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( LDCCaloDigi_srcs ${LDCCaloDigi_cxx_srcs} )

####################################################################################################
# ClusterCheater

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/ClusterCheater/src  ClusterCheater_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ClusterCheater_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( ClusterCheater_srcs ${ClusterCheater_cxx_srcs} )

####################################################################################################
# ClusterCheater5_3

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/ClusterCheater5_3/src  ClusterCheater5_3_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ClusterCheater5_3_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( ClusterCheater5_3_srcs ${ClusterCheater5_3_cxx_srcs} )

####################################################################################################
# NNClustering

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/NNClustering/src NNClustering_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${NNClustering_cxx_srcs} PROPERTIES
    # require proper c++
    # FIXME -pedantic flag crashes under sl4
    COMPILE_FLAGS "-Wall -ansi" )

SET( NNClustering_srcs ${NNClustering_cxx_srcs} )

####################################################################################################
# TrackwiseClustering

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/TrackwiseClustering/src TrackwiseClustering_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackwiseClustering_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( TrackwiseClustering_srcs ${TrackwiseClustering_cxx_srcs} )

####################################################################################################
# PhotonFinderKit

# cxx sources
AUX_SOURCE_DIRECTORY( ../Clustering/PhotonFinderKit/src PhotonFinderKit_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${PhotonFinderKit_cxx_srcs} PROPERTIES
    # require proper c++
    # FIXME -pedantic flag crashes under sl4
    COMPILE_FLAGS "-Wall -ansi" )

SET( PhotonFinderKit_srcs ${PhotonFinderKit_cxx_srcs} )

####################################################################################################
# Wolf

# cxx sources
AUX_SOURCE_DIRECTORY( ../Pflow/Wolf/src Wolf_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${Wolf_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( Wolf_srcs ${Wolf_cxx_srcs} )

####################################################################################################
# TrackBasedPFlow

# cxx sources
AUX_SOURCE_DIRECTORY( ../Pflow/TrackBasedPFlow/src TrackBasedPFlow_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackBasedPFlow_cxx_srcs} PROPERTIES
    # require proper c++
    # FIXME -pedantic flag crashes under sl4
    COMPILE_FLAGS "-Wall -ansi" )

SET( TrackBasedPFlow_srcs ${TrackBasedPFlow_cxx_srcs} )

####################################################################################################
# FTDDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/FTDDigi/src FTDDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${FTDDigi_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( FTDDigi_srcs ${FTDDigi_cxx_srcs} )

####################################################################################################
# ETDDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/ETDDigi/src ETDDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${ETDDigi_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( ETDDigi_srcs ${ETDDigi_cxx_srcs} )

####################################################################################################
# TPCDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/TPCDigi/src TPCDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TPCDigi_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( TPCDigi_srcs ${TPCDigi_cxx_srcs} )

####################################################################################################
# VTXDigi

# cxx sources
AUX_SOURCE_DIRECTORY( ../TrackDigi/VTXDigi/src VTXDigi_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${VTXDigi_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( VTXDigi_srcs ${VTXDigi_cxx_srcs} )

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
    COMPILE_FLAGS "-ansi -pedantic" )

# fortran sources
FILE( GLOB BrahmsTracking_f77_srcs ../Tracking/BrahmsTracking/src/f77/*.F )

SET( BrahmsTracking_srcs ${BrahmsTracking_cxx_srcs} ${BrahmsTracking_c_srcs} ${BrahmsTracking_f77_srcs} )

####################################################################################################
# FullLDCTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/FullLDCTracking/src FullLDCTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${FullLDCTracking_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( FullLDCTracking_srcs ${FullLDCTracking_cxx_srcs} )

####################################################################################################
# SiliconTracking

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/SiliconTracking/src SiliconTracking_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${SiliconTracking_cxx_srcs} PROPERTIES
    # require proper c++
    # FIXME -pedantic
    COMPILE_FLAGS "-Wall -ansi" )

SET( SiliconTracking_srcs ${SiliconTracking_cxx_srcs} )

####################################################################################################
# TrackCheater

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/TrackCheater/src TrackCheater_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${TrackCheater_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( TrackCheater_srcs ${TrackCheater_cxx_srcs} )

####################################################################################################
# V0Finder

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/V0Finder/src V0Finder_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${V0Finder_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( V0Finder_srcs ${V0Finder_cxx_srcs} )

####################################################################################################
# KinkFinder

# cxx sources
AUX_SOURCE_DIRECTORY( ../Tracking/KinkFinder/src KinkFinder_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${KinkFinder_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( KinkFinder_srcs ${KinkFinder_cxx_srcs} )


####################################################################################################
# PFOID

# cxx sources
AUX_SOURCE_DIRECTORY( ../PFOID/src PFOID_cxx_srcs )

SET_SOURCE_FILES_PROPERTIES( ${PFOID_cxx_srcs} PROPERTIES
    # require proper c++
    COMPILE_FLAGS "-Wall -ansi -pedantic" )

SET( PFOID_srcs ${PFOID_cxx_srcs} )


####################################################################################################
### LIBRARY
####################################################################################################

# set library sources
SET( lib_srcs )
FOREACH( pkgDir ${MarlinReco_PKGS} )
    # extract basename from dir (e.g. Analysis/EventShapes -> EventShapes)
    STRING( REGEX REPLACE ".*/([^ ]+)" "\\1" pkgName ${pkgDir} )
    # all c++ sources need -Wno-long-long flag (due to lcio)
    SET_SOURCE_FILES_PROPERTIES( ${${pkgName}_cxx_srcs} PROPERTIES COMPILE_FLAGS "-Wno-long-long" )
    SET( lib_srcs ${lib_srcs} ${${pkgName}_srcs} )
ENDFOREACH()

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# workaround for fortran cmake bug in cmake 2.4.x versions
IF( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.6 )
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
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > TE.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > helixextrapolation.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > old.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > bit.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > copy.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > polar.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > primary.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > theta.mod.proxy
    
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Fixing fortran bug for cmake 2.4 versions..." VERBATIM )
    
    # establish dependency for fixing fortran bug
    ADD_DEPENDENCIES( lib_${PROJECT_NAME} fix_f77_cmake_bug )
ENDIF()
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

SET( G2C_LIB G2C_LIB-NOTFOUND )
# prevent static library to be detected instead of shared one
FIND_LIBRARY( G2C_LIB NAMES ${G2C_LIBNAME} PATHS ${G2C_LIB_PATHS} NO_DEFAULT_PATH )
# if previous find fails, try to find in system paths..
FIND_LIBRARY( G2C_LIB NAMES ${G2C_LIBNAME} )

IF( NOT G2C_LIB )
    MESSAGE( FATAL_ERROR "Check for lib${G2C_LIBNAME}: failed to find ${G2C_LIBNAME} library!!" )
ELSE()
    MESSAGE( STATUS "Check for lib${G2C_LIBNAME}: ${G2C_LIB}" )
ENDIF()
MARK_AS_ADVANCED( G2C_LIB )

# link MarlinReco library to libg2c/libgfortran
TARGET_LINK_LIBRARIES( lib_${PROJECT_NAME} ${G2C_LIB} )
IF(APPLE)
  TARGET_LINK_LIBRARIES( lib_${PROJECT_NAME} streamlog )
ENDIF(APPLE)

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )

