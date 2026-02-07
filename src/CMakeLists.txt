target_sources(${PROJECT_NAME} PRIVATE
  AnalysisWaveform.cc
  BandwidthMeasure.cc
  BasicFilters.cc
  CorrelationSummary.cc
  DiodeFilter.cc
  EventFaker.cc
  EventSummary.cc
  FilterOperation.cc
  FilterStrategy.cc
  FilteredEvent.cc
  FreqDomainFunction.cc
  GeomFilter.cc
  ImpulsivityMeasure.cc
  InterpolatedNotch.cc
  NoiseMonitor.cc
  NoiseSummary.cc
  Polarimetry.cc
  PrettyAnalysisWaveform.cc
  ResponseManager.cc
  SensitivityCalculator.cc
  SystemResponse.cc
  TGraphAligned.cc
  TMVA.cc
  TaperFilter.cc
  Templates.cc
  UsefulAttitude.cc
)
set(HEADER_FILES
  pueo/AnalysisWaveform.h
  pueo/BandwidthMeasure.h
  pueo/BasicFilters.h
  pueo/CorrelationSummary.h
  pueo/DiodeFilter.h
  pueo/EventFaker.h
  pueo/EventReconstructor.h
  pueo/EventSummary.h
  pueo/FilterOperation.h
  pueo/FilterStrategy.h
  pueo/FilteredEvent.h
  pueo/FreqDomainFunction.h
  pueo/GeomFilter.h
  pueo/ImpulsivityMeasure.h
  pueo/InterpolatedNotch.h
  pueo/NoiseMonitor.h
  pueo/NoiseSummary.h
  pueo/Polarimetry.h
  pueo/PrettyAnalysisWaveform.h
  pueo/ResponseManager.h
  pueo/SensitivityCalculator.h
  pueo/SystemResponse.h
  pueo/TGraphAligned.h
  pueo/TMVA.h
  pueo/TaperFilter.h
  pueo/Templates.h
  pueo/UsefulAttitude.h
)

target_sources(${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  BASE_DIRS ${CMAKE_CURRENT_LIST_DIR}
  FILES ${HEADER_FILES}
)

# CERN ROOT dictionary generation
set(DICTNAME G__${PROJECT_NAME})
root_generate_dictionary(${DICTNAME} ${HEADER_FILES} MODULE ${PROJECT_NAME} LINKDEF LinkDef.h)

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

