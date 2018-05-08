get_filename_component(TARGET ${CMAKE_CURRENT_SOURCE_DIR} NAME)
FILE(GLOB SRC_FILES *.cpp *.h)
SET(TARGET_SRC ${SRC_FILES})

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.5.0 QUIET HINTS ${deal.II_DIR} ${DEAL_II_DIR})
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
include_directories(../..)
add_executable(${TARGET} ${TARGET_SRC})
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
DEAL_II_SETUP_TARGET(${TARGET} DEBUG)
ELSE()
DEAL_II_SETUP_TARGET(${TARGET} RELEASE)
ENDIF()
target_link_libraries(${TARGET} mhdeal)