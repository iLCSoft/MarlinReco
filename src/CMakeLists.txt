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

        # BCalTagEfficiency fix
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > sine.mod.proxy
        COMMAND ${CMAKE_COMMAND} -E echo "dummy file: workaround for fortran cmake bug" > "...).mod.proxy"
    
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

IF( NOT BUILD_32BIT_COMPATIBLE )
    SET( FORTRAN_LIB_PATHS
        # 64 bit paths
        /usr/lib/gcc/x86_64-redhat-linux/3.4.3      # SL4 native 64bit
        /usr/lib/gcc/x86_64-redhat-linux/3.4.6      # SL5 native 64bit
        /usr/lib/gcc/x86_64-redhat-linux/4.1.1      # SL5 native 64bit
        /usr/lib64
    )
ENDIF()

SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        # 32 bit paths
        /usr/lib/gcc-lib/i386-redhat-linux/3.2.3    # SL3 native 32bit
        /usr/lib/gcc/i386-redhat-linux/3.4.3        # SL4 native 32bit
        /usr/lib/gcc/i386-redhat-linux/3.4.6        # SL5 native 32bit
)

IF( BUILD_32BIT_COMPATIBLE )
    SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        /usr/lib/gcc/x86_64-redhat-linux/3.4.3/32   # SL4 64bit 32bit compatible
        /usr/lib/gcc/x86_64-redhat-linux/3.4.6/32   # SL5 64bit 32bit compatible
        /usr/lib/gcc/x86_64-redhat-linux/4.1.1/32   # SL5 64bit 32bit compatible
        /usr/lib
    )
ENDIF()

SET( FORTRAN_LIB_PATHS ${FORTRAN_LIB_PATHS}
        /usr/lib/gcc/x86_64-linux-gnu/4.4           # Ubuntu 10.04
        /usr/osxws/lib                              # OSX
)

MESSAGE( STATUS "FORTRAN_LIB_PATHS: ${FORTRAN_LIB_PATHS}" )

SET( FORTRAN_LIB FORTRAN_LIB-NOTFOUND )
FIND_LIBRARY( FORTRAN_LIB NAMES ${FORTRAN_LIBNAME} PATHS ${FORTRAN_LIB_PATHS} NO_DEFAULT_PATH )
# if previous find fails, try to find in system paths..
FIND_LIBRARY( FORTRAN_LIB NAMES ${FORTRAN_LIBNAME} )

IF( NOT FORTRAN_LIB )
    MESSAGE( FATAL_ERROR "failed to find ${FORTRAN_LIBNAME} library!!" )
ELSE()
    MESSAGE( STATUS "Check for lib${FORTRAN_LIBNAME}: ${FORTRAN_LIB}" )
ENDIF()
MARK_AS_ADVANCED( FORTRAN_LIB )

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

