SET(TARGET "mhdeal")

SET(TARGET_SRC
  util.h
  equations.h
  equationsMhd.cpp
  DealiiExtensions.cpp
  DealiiExtensions.h
  fe_div_free.h
  fe_div_free.cpp
  fe_taylor.h
  fe_taylor.cpp
  fe_is_constant_interface.h
  equationsMhd.h
  parameters.h
  boundaryConditions.cpp
  boundaryConditions.h
  initialCondition.cpp
  initialCondition.h
  problem.cpp
  problem.h
  complete_elliptic_integrals.cc
  complete_elliptic_integrals.h
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
include_directories(${DEAL_II_INCLUDE_DIRS})
add_library(${TARGET} ${TARGET_SRC})
DEAL_II_SETUP_TARGET(${TARGET} RELEASE)
add_subdirectory(examples)
install(TARGETS ${TARGET} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)