cmake_minimum_required(VERSION 3.23.0)
project(pueoAnalysisFramework VERSION 0.0.1)
add_library(${PROJECT_NAME} SHARED)

target_sources(${PROJECT_NAME} PRIVATE
  # src/AnalysisWaveform.cc
  # src/BandwidthMeasure.cc
  # src/BasicFilters.cc
  # src/CorrelationSummary.cc
  # src/DiodeFilter.cc
  # src/EventFaker.cc
  # src/EventSummary.cc
  # src/FilterOperation.cc
  # src/FilterStrategy.cc
  # src/FilteredEvent.cc
  # src/FreqDomainFunction.cc
  # src/GeomFilter.cc
  # src/ImpulsivityMeasure.cc
  # src/InterpolatedNotch.cc
  # src/NoiseMonitor.cc
  # src/NoiseSummary.cc
  # src/Polarimetry.cc
  # src/PrettyAnalysisWaveform.cc
  # src/ResponseManager.cc
  # src/SensitivityCalculator.cc
  # src/SystemResponse.cc
  # src/TGraphAligned.cc
  # src/TMVA.cc
  # src/TaperFilter.cc
  # src/Templates.cc
  # src/UsefulAttitude.cc
)
set(HEADER_FILES
  # src/pueo/AnalysisWaveform.h
  # src/pueo/BandwidthMeasure.h
  # src/pueo/BasicFilters.h
  # src/pueo/CorrelationSummary.h
  # src/pueo/DiodeFilter.h
  # src/pueo/EventFaker.h
  # src/pueo/EventReconstructor.h
  src/pueo/EventSummary.h
  # src/pueo/FilterOperation.h
  # src/pueo/FilterStrategy.h
  # src/pueo/FilteredEvent.h
  # src/pueo/FreqDomainFunction.h
  # src/pueo/GeomFilter.h
  # src/pueo/ImpulsivityMeasure.h
  # src/pueo/InterpolatedNotch.h
  # src/pueo/NoiseMonitor.h
  # src/pueo/NoiseSummary.h
  # src/pueo/Polarimetry.h
  # src/pueo/PrettyAnalysisWaveform.h
  # src/pueo/ResponseManager.h
  # src/pueo/SensitivityCalculator.h
  # src/pueo/SystemResponse.h
  # src/pueo/TGraphAligned.h
  # src/pueo/TMVA.h
  # src/pueo/TaperFilter.h
  # src/pueo/Templates.h
  # src/pueo/UsefulAttitude.h
)

target_sources(${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  BASE_DIRS src  # <-- Everything in BASE_DIRS is available during the build
  FILES ${HEADER_FILES}  # <-- But only FILES will be installed, although these two are currently the same
)


## ------------- Dependency -------------

# find_package(ROOT CONFIG REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net TMVA RooFit RooStats) 
find_package(ROOT CONFIG REQUIRED COMPONENTS TreePlayer)
find_package(AntarcticaRoot CONFIG REQUIRED) # this "transitively" fetches RootFftwWrapper
find_package(pueoEvent CONFIG REQUIRED)
find_package(OpenMP)
## todo, OpenMP can probably also come from RootFftwWrapper but not sure how much that simplifies things.
#  in the meantime, references:
#    1) https://cmake.org/cmake/help/latest/module/FindOpenMP.html#result-variables
#    2) https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PRIVATE USE_OMP)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
else()
  message(STATUS "OpenMP not found, compiling without it.")
endif()

target_link_libraries(${PROJECT_NAME} 
  PUBLIC ROOT::TreePlayer AntarcticaRoot::AntarcticaRoot PUEO::pueoEvent
)
# CERN ROOT dictionary generation and installation
set(DICTNAME G__${PROJECT_NAME})
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION lib
)



## ------------- Installation -------------

include(GNUInstallDirs)
# the config file has instruction for downstream projects on how to fetch dependencies of this project
install(
  FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)

install(
  TARGETS ${PROJECT_NAME} # just one target, ie this library, for now
  EXPORT ${PROJECT_NAME}Targets # put all targets listed above with this particular "export"
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # install the header files of the target
  LIBRARY DESTINATION lib
)

install(
  EXPORT ${PROJECT_NAME}Targets # auto create and install pueoAnalysisFrameworkTargets.cmake directly to install tree
  NAMESPACE PUEO::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}
)



## ------------- Options -------------

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelwithDebugInfo CACHE STRING "" FORCE)
  # use warning flags in RelWithDebInfo mode, which is set to be the default mode
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:RelWithDebInfo>:-Wall -Wextra -Wno-unused-function>)
endif()

## Change Log Feb 8, 2026
#  The option "VECTORIZE" is removed.
#  The relevant additional include path and linking should come from the upstream package RootFftwWrapper.
#  There, the option is on by default (although pueoBuilder turns this off automatically for arm64 platforms).
#  That is, the compiler flag "-DENABLE_VECTORIZE" is public and will propagate from RootFftwWrapper 
#  to this project; the include path and linkage also comes from there.


## ------------- TODO -------------

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


## todo, add this back if it's still needed
# option( PUEO_ANALYSIS_DEBUG "Enable some debug options that probably slow things down " )
# if (PUEO_ANALYSIS_DEBUG)
#   add_definitions( -DPUEO_ANALYSIS_DEBUG)
# endif()
