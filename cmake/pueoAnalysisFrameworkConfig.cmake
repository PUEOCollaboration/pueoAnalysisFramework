include("${CMAKE_CURRENT_LIST_DIR}/pueoAnalysisFrameworkTargets.cmake")

include(CMakeFindDependencyMacro)
find_dependency(AntarcticaRoot CONFIG REQUIRED)
find_dependency(pueoEvent CONFIG REQUIRED)
find_dependency(ROOT CONFIG REQUIRED COMPONENTS TreePlayer)

