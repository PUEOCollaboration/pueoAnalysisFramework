
# CMakeLists.txt for pueoAnalysisFramework
cmake_minimum_required(VERSION 3.0.0)
project(pueoAnalysisFramework)

## perhaps this should be changed?
set(libname "pueoAnalysis")


# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


if(DEFINED ENV{PUEO_UTIL_INSTALL_DIR})
    message("PUEO_UTIL_INSTALL_DIR is set to $ENV{PUEO_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{PUEO_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{PUEO_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  PUEO_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UTIL_SHARE_DIR}/cmake/modules/")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
option(ENABLE_ROOSTATS "Enable RooStats/RooFit") 
if (ENABLE_ROOSTATS) 
  find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net TMVA RooFit RooStats) 
  add_definitions(-DUSE_ROOSTATS) 
  message("Using RooStats") 
else() 
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net TMVA)
endif() 

find_package(FFTW REQUIRED)

link_directories(${UTIL_LIB_DIR})


include(${ROOT_USE_FILE})

message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")


include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})

file(GLOB HEADER_FILES
    "src/pueo/*.h"
    )

file(GLOB_RECURSE SOURCE_FILES src/*.cc)

set(DICTNAME G__${libname})

#-- Wavelet Deconvolution using FoRWarD.
option(WAVELETDECONV "Enable Wavelet Deconvolution." OFF)
if (WAVELETDECONV)

  # the directory where we store the FoRWarD wrapper
  set(FORWARD_DIR "${PROJECT_SOURCE_DIR}/thirdparty/forward/anita/")

  # due to a cylic dependency, we inject forwaord directly
  # into the source and header tooling
  list(APPEND SOURCE_FILES "${FORWARD_DIR}/WaveletDeconvolution.cc")
  list(APPEND HEADER_FILES "${FORWARD_DIR}/WaveletDeconvolution.h")

  # and then use the forward directory as an include
  include_directories("${PROJECT_SOURCE_DIR}/thirdparty/forward/include")

  # and then use the forward directory as an include
  include_directories("${FORWARD_DIR}")

  # and install FoRWarD into PUEO_UTIL_INSTALL_DIR
  install(DIRECTORY
    "${PROJECT_SOURCE_DIR}/thirdparty/forward/include/forward"
    DESTINATION ${UTIL_INC_DIR})

endif()

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with generated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ZLIB_LIBRARIES} RootFftwWrapper pueoEvent AntarcticaRoot ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})


message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")
message("DICTNAME is set to ${DICTNAME}" )

install (FILES ${HEADER_FILES} DESTINATION ${UTIL_INC_DIR}/pueo)
install (TARGETS ${libname}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
  install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
#  install (FILES ${PROJECT_BINARY_DIR}/lib${libname}.rootmap DESTINATION ${UTIL_LIB_DIR})
endif()

#################################################################################3
### The default build configuration is INSANE. No optimization? What is this 1971?

message (" Adding new build type")

set(CMAKE_CXX_FLAGS_DEFAULT
  "-Os -g -Wall"
  CACHE STRING "c++ Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_C_FLAGS_DEFAULT
  "-Os -g -Wall"
  CACHE STRING "c Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_EXE_LINKER_FLAGS_DEFAULT
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds"
  FORCE )

set(CMAKE_SHARED_LINKER_FLAGS_DEFAULT
  "-g"
  CACHE STRING "ld Flags used during default eventReaderRoot builds"
  FORCE )


mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose tye type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE )
endif()

#################################################################################
###  Manual vectorization option
###  Maybe in the future this can be made a CMake module or something
#################################################################################

option(VECTORIZE "Enable Manual SIMD Vectorization" ON) 
set (VECTORDIR ${UTIL_INC_DIR}/vectorclass)
if (VECTORIZE)
  add_definitions( -DENABLE_VECTORIZE )

  if(CMAKE_COMPILER_IS_GNUCXX)

    ### someone should do this for clang if they want it to be as fast as possible
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.1)
      if(NATIVE_ARCH) #define in anitaBuildTool, default on
	add_definitions( -march=native)
      endif()
    endif()
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
      # Vectorize docs...
      # If you are using the Gnu compiler version 3.x or 4.x then you must
      # set the ABI version to 4 or more, or 0 for a reasonable default
      add_definitions(-fabi-version=0)
    endif()
  endif()

  #  message(${VECTORDIR})
  include_directories(${VECTORDIR}) 


endif (VECTORIZE)

option ( FORCE_OLD_GPP_ABI " Force old g++ ABI; this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations" OFF )

if (FORCE_OLD_GPP_ABI)
  add_definitions( -D_GLIBCXX_USE_CXX11_ABI=0 )
endif()

option( ENABLE_OPENMP "Enable OpenMP support")



if(ENABLE_OPENMP)
  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    add_definitions (-DUSE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support")
  endif()
endif()

install (DIRECTORY data/ DESTINATION ${UTIL_SHARE_DIR}/pueoAnalysisFramework)

option( PUEO_ANALYSIS_DEBUG "Enable some debug options that probably slow things down " )


if (PUEO_ANALYSIS_DEBUG)
  add_definitions( -DPUEO_ANALYSIS_DEBUG)
endif()
