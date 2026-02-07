cmake_minimum_required(VERSION 3.23.0)
project(pueoAnalysisFramework VERSION 0.0.1)
add_library(${PROJECT_NAME} SHARED)

# find_package(ROOT CONFIG REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net TMVA RooFit RooStats) 
find_package(ROOT CONFIG REQUIRED COMPONENTS TreePlayer)
find_package(AntarcticaRoot CONFIG REQUIRED) # this "transitively" fetches RootFftwWrapper
find_package(pueoEvent CONFIG REQUIRED)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelwithDebugInfo CACHE STRING "" FORCE)
  # use warning flags in RelWithDebInfo mode, which is set to be the default mode
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:RelWithDebInfo>:-Wall -Wextra -Wno-unused-function>)
endif()

add_subdirectory(src) # this calls src/CMakeLists.txt to add cpp and header files

target_link_libraries(${PROJECT_NAME} 
  PUBLIC ROOT::TreePlayer AntarcticaRoot::AntarcticaRoot PUEO::pueoEvent
)

# the config file has instruction for downstream projects on how to fetch dependencies of this project
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

install(
  TARGETS ${PROJECT_NAME} # just one target, ie this library, for now
  EXPORT ${PROJECT_NAME}Targets # put all targets listed above with this particular "export"
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # install the header files of the target
)

install(
  EXPORT ${PROJECT_NAME}Targets # auto create and install pueoAnalysisFrameworkTargets.cmake directly to install tree
  NAMESPACE PUEO::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)


## todo
# #-- Wavelet Deconvolution using FoRWarD.
# option(WAVELETDECONV "Enable Wavelet Deconvolution." OFF)
# if (WAVELETDECONV)
#
#   # the directory where we store the FoRWarD wrapper
#   set(FORWARD_DIR "${PROJECT_SOURCE_DIR}/thirdparty/forward/anita/")
#
#   # due to a cylic dependency, we inject forwaord directly
#   # into the source and header tooling
#   list(APPEND SOURCE_FILES "${FORWARD_DIR}/WaveletDeconvolution.cc")
#   list(APPEND HEADER_FILES "${FORWARD_DIR}/WaveletDeconvolution.h")
#
#   # and then use the forward directory as an include
#   include_directories("${PROJECT_SOURCE_DIR}/thirdparty/forward/include")
#
#   # and then use the forward directory as an include
#   include_directories("${FORWARD_DIR}")
#
#   # and install FoRWarD into PUEO_UTIL_INSTALL_DIR
#   install(DIRECTORY
#     "${PROJECT_SOURCE_DIR}/thirdparty/forward/include/forward"
#     DESTINATION ${UTIL_INC_DIR})
#
# endif()


## todo, this can come from RootFftwWrapper
# option(VECTORIZE "Enable Manual SIMD Vectorization" ON) 
# set (VECTORDIR ${UTIL_INC_DIR}/vectorclass)
# if (VECTORIZE)
#   add_definitions( -DENABLE_VECTORIZE )
#
#   if(CMAKE_COMPILER_IS_GNUCXX)
#
#     ### someone should do this for clang if they want it to be as fast as possible
#     if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.1)
#       if(NATIVE_ARCH) #define in anitaBuildTool, default on
# 	add_definitions( -march=native)
#       endif()
#     endif()
#     if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
#       # Vectorize docs...
#       # If you are using the Gnu compiler version 3.x or 4.x then you must
#       # set the ABI version to 4 or more, or 0 for a reasonable default
#       add_definitions(-fabi-version=0)
#     endif()
#   endif()
#
#   #  message(${VECTORDIR})
#   include_directories(${VECTORDIR}) 
# endif (VECTORIZE)


## todo, this can once again come from RootFftwWrapper
# option( ENABLE_OPENMP "Enable OpenMP support")
# if(ENABLE_OPENMP)
#   FIND_PACKAGE( OpenMP REQUIRED)
#   if(OPENMP_FOUND)
#     message("OPENMP FOUND")
#     add_definitions (-DUSE_OMP)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#   else()
#     message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support")
#   endif()
# endif()

## todo, add this back if it's still needed
# option( PUEO_ANALYSIS_DEBUG "Enable some debug options that probably slow things down " )
# if (PUEO_ANALYSIS_DEBUG)
#   add_definitions( -DPUEO_ANALYSIS_DEBUG)
# endif()
