########################################################
# cmake file for building MarlinReco documentation
# @author Jan Engels, DESY
########################################################


FIND_PACKAGE( Doxygen )
FIND_PACKAGE( LATEX )

SET( DOC_SRC_DIR "${PROJECT_SOURCE_DIR}/doc" )
SET( DOC_BIN_DIR "${PROJECT_BINARY_DIR}/docbuild" )
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${DOC_BIN_DIR}" )

# global doc target for calling doc_html and doc_manual
ADD_CUSTOM_TARGET( doc )

# doxygen
IF( DOXYGEN_FOUND )

    #SET( DOX_INPUT overview.html ../include ../src )
    SET( DOX_PREDEFINED MARLIN_USE_AIDA ) # TODO test !!
    
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_BIN_DIR}/html/index.html"
        COMMAND DOX_PROJECT_NAME="${PROJECT_NAME}"
                DOX_PROJECT_NUMBER="${${PROJECT_NAME}_VERSION}"
                DOX_OUTPUT_DIRECTORY="${DOC_BIN_DIR}"
                #DOX_INPUT="${DOX_INPUT}"
                DOX_PREDEFINED="${DOX_PREDEFINED}"
                "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${DOC_SRC_DIR}"
        COMMENT "Building API Documentation..."
        DEPENDS Doxyfile CMakeLists.txt overview.html
    )
    
    # add doc_html target
    ADD_CUSTOM_TARGET( doc_html DEPENDS "${DOC_BIN_DIR}/html/index.html" )

    # tell doc target to call doc_html
    ADD_DEPENDENCIES( doc doc_html )
    
ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC set to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()



# MANUAL
FIND_PACKAGE( LATEX )
IF( LATEX_COMPILER AND MAKEINDEX_COMPILER AND PDFLATEX_COMPILER
    AND DVIPS_CONVERTER AND LATEX2HTML_CONVERTER )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_BIN_DIR}/manual/manual.tex"
        COMMAND "${CMAKE_COMMAND}" -E copy_directory "${DOC_SRC_DIR}/manual" "${DOC_BIN_DIR}/manual"
        DEPENDS ${DOC_SRC_DIR}/manual/*
    )

    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_BIN_DIR}/manual_html/manual.html"
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${LATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        #COMMAND "${PDFLATEX_COMPILER}" ARGS -interaction=batchmode manual.tex
        COMMAND "${LATEX2HTML_CONVERTER}" ARGS -mkdir -dir ../manual_html -verbosity 0 -info 0 -no_auto_link -split 0 -no_navigation manual.tex
        COMMAND "${DVIPS_CONVERTER}" ARGS -q -o manual.ps manual.dvi
        WORKING_DIRECTORY "${DOC_BIN_DIR}/manual"
        COMMENT "Building Documentation - Users Manual..."
        DEPENDS "${DOC_BIN_DIR}/manual/manual.tex"
    )

    # add doc_manual target
    ADD_CUSTOM_TARGET( doc_manual DEPENDS "${DOC_BIN_DIR}/manual_html/manual.html" )

    # tell doc target to call doc_manual
    ADD_DEPENDENCIES( doc doc_manual )
    
ELSE()
    MESSAGE( STATUS "Latex not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC set to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# install documentation
IF( INSTALL_DOC )

    # execute 'make doc' before 'make install'
    INSTALL( CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_BUILD_TOOL} doc)" )

    # in-source vs out-of-source installations
    IF( CMAKE_INSTALL_PREFIX STREQUAL "${PROJECT_SOURCE_DIR}" )

        # --- in-source installations ---
        SET( DOC_INSTALL_PREFIX doc )

        ## keep source tree clean
        ADD_CUSTOM_TARGET( uninstall_doc
            COMMAND "${CMAKE_COMMAND}" -E remove_directory "${DOC_SRC_DIR}/html"
            COMMAND "${CMAKE_COMMAND}" -E remove_directory "${DOC_SRC_DIR}/manual_html"
        )
        ADD_DEPENDENCIES( uninstall uninstall_doc )

    ELSE()

        # --- out-of-source installations ---
        SET( DOC_INSTALL_PREFIX doc/${PROJECT_NAME} )

        INSTALL( FILES release.notes FullLDCTracking_Manual.pdf DESTINATION ${DOC_INSTALL_PREFIX} )

    ENDIF( CMAKE_INSTALL_PREFIX STREQUAL "${PROJECT_SOURCE_DIR}" )

    # install html documentation
    INSTALL_DIRECTORY( "${DOC_BIN_DIR}/html" DESTINATION ${DOC_INSTALL_PREFIX} )
    INSTALL_DIRECTORY( "${DOC_BIN_DIR}/manual_html" DESTINATION ${DOC_INSTALL_PREFIX} )

ENDIF()

