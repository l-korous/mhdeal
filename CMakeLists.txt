SET(TARGET "mhdeal")

SET(TARGET_SRC
  adaptivity.cpp
  adaptivity.h
  assemblingUtilities.cpp
  assemblingUtilities.h
  util.h
  equations.h
  equationsMhd.cpp
  feIsConstInterface.h
  dealiiExtensions.cpp
  dealiiExtensions.h
  feDivFree.h
  feDivFree.cpp
  feTaylor.h
  feTaylor.cpp
  equationsMhd.h
  parameters.cpp
  parameters.h
  boundaryConditions.cpp
  boundaryConditions.h
  initialCondition.cpp
  initialCondition.h
  problem.cpp
  problem.h
  completeEllipticIntegrals.cpp
  completeEllipticIntegrals.h
  numericalFlux.cpp
  numericalFlux.h
  slopeLimiter.cpp
  slopeLimiter.h
)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Are all dependencies fulfilled?
#
IF(NOT DEAL_II_WITH_TRILINOS)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_TRILINOS = ON
One or all of these are OFF in your installation but are required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
IF(WIN32)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LTCG")
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Od /MP")
    ELSE()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ob2 /Oi /Ot /GL /MP")
    ENDIF()
ELSE()
    IF(NOT(CMAKE_BUILD_TYPE STREQUAL "Debug"))
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    ENDIF()
ENDIF()
include_directories(${DEAL_II_INCLUDE_DIRS})
add_library(${TARGET} ${TARGET_SRC})
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
DEAL_II_SETUP_TARGET(${TARGET} DEBUG)
ELSE()
DEAL_II_SETUP_TARGET(${TARGET} RELEASE)
ENDIF()
add_subdirectory(examples)
install(TARGETS ${TARGET} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)