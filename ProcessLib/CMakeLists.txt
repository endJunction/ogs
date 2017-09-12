include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES CalculateSurfaceFlux)
APPEND_SOURCE_FILES(SOURCES Deformation)
APPEND_SOURCE_FILES(SOURCES Parameter)
APPEND_SOURCE_FILES(SOURCES PhaseField)
APPEND_SOURCE_FILES(SOURCES SmallDeformation)
APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    PUBLIC BaseLib GeoLib MaterialLib MathLib MeshLib NumLib logog
    PRIVATE MeshGeoToolsLib
)

ADD_VTK_DEPENDENCY(ProcessLib)

if(TARGET Eigen)
    add_dependencies(ProcessLib Eigen)
endif()
if(TARGET Boost)
    add_dependencies(ProcessLib Boost)
endif()
if(OGS_INSITU)
    target_link_libraries(ProcessLib InSituLib)
endif()

if(OGS_USE_PCH)
    cotire(ProcessLib)
endif()

### Tests ###
include(ComponentTransport/Tests.cmake)
include(GroundwaterFlow/Tests.cmake)
include(HeatConduction/Tests.cmake)
include(HT/Tests.cmake)
include(HydroMechanics/Tests.cmake)
include(LIE/HydroMechanics/Tests.cmake)
include(LIE/SmallDeformation/Tests.cmake)
include(LiquidFlow/Tests.cmake)
include(PhaseField/Tests.cmake)
include(RichardsFlow/Tests.cmake)
include(SmallDeformation/Tests.cmake)
include(TES/Tests.cmake)
include(ThermoMechanics/Tests.cmake)
include(TwoPhaseFlowWithPP/Tests.cmake)
include(TwoPhaseFlowWithPrho/Tests.cmake)
include(ThermalTwoPhaseFlowWithPP/Tests.cmake)
