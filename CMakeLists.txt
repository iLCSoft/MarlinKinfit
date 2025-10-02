########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 3.23 FATAL_ERROR)
CMAKE_POLICY(SET CMP0033 OLD)  # Allow to call export_library_dependencies
########################################################


# project name
PROJECT( MarlinKinfit )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 6 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )

SET( ${PROJECT_NAME}_VERSION  "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}" )


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( LCIO 2.21 REQUIRED )
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( GSL 1.12 REQUIRED )
FIND_PACKAGE( ROOT 6.14 REQUIRED )

FIND_PACKAGE( AIDA )
IF( AIDA_FOUND )
    MESSAGE( STATUS "AIDA -- found" )
ELSE()
    MESSAGE( STATUS "AIDA -- not found" )
ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################


# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
target_link_libraries( ${PROJECT_NAME} PUBLIC
    LCIO::lcio
    GSL::gsl
    CLHEP::CLHEP
    ROOT::Core ROOT::Physics ROOT::Hist ROOT::Tree
)
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
)
if (AIDA_FOUND)
    target_link_libraries( ${PROJECT_NAME} PUBLIC
        ${AIDA_LIBARIES}
    )
    target_compile_definitions( ${PROJECT_NAME} PUBLIC -DMARLIN_USE_AIDA)
    target_include_directories( ${PROJECT_NAME} PUBLIC ${AIDA_INCLUDE_DIRS} )
endif()


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# Install library
# include directories
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT shlib
)

include(cmake/createMarlinKinfitConfig.cmake)
