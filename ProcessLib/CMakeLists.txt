include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES CalculateSurfaceFlux)
APPEND_SOURCE_FILES(SOURCES ComponentTransport)
APPEND_SOURCE_FILES(SOURCES Deformation)
APPEND_SOURCE_FILES(SOURCES GroundwaterFlow)
APPEND_SOURCE_FILES(SOURCES HeatConduction)
APPEND_SOURCE_FILES(SOURCES HT)
APPEND_SOURCE_FILES(SOURCES HydroMechanics)
if(Protobuf_FOUND)
endif() # Protobuf_FOUND
APPEND_SOURCE_FILES(SOURCES LIE/BoundaryCondition)
APPEND_SOURCE_FILES(SOURCES LIE/Common)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics)
APPEND_SOURCE_FILES(SOURCES LIE/HydroMechanics/LocalAssembler)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation)
APPEND_SOURCE_FILES(SOURCES LIE/SmallDeformation/LocalAssembler)
APPEND_SOURCE_FILES(SOURCES LiquidFlow)
APPEND_SOURCE_FILES(SOURCES Parameter)
APPEND_SOURCE_FILES(SOURCES RichardsFlow)

if(Protobuf_FOUND)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS SmallDeformationCommon/integration_point.proto)
list(APPEND SOURCES ${PROTO_SRCS} ${PROTO_HDRS})
endif() # Protobuf_FOUND

APPEND_SOURCE_FILES(SOURCES SmallDeformation)
APPEND_SOURCE_FILES(SOURCES SmallDeformationNonlocal)
APPEND_SOURCE_FILES(SOURCES TES)
APPEND_SOURCE_FILES(SOURCES TwoPhaseFlowWithPP)
APPEND_SOURCE_FILES(SOURCES TwoPhaseFlowWithPrho)
APPEND_SOURCE_FILES(SOURCES Utils)

add_library(ProcessLib ${SOURCES})

target_link_libraries(ProcessLib
    PUBLIC BaseLib GeoLib MaterialLib MathLib MeshLib NumLib logog
    PRIVATE MeshGeoToolsLib
)

if(Protobuf_FOUND)
    target_link_libraries(ProcessLib ${PROTOBUF_LIBRARIES})
endif() # Protobuf_FOUND

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
include(RichardsFlow/Tests.cmake)
include(SmallDeformation/Tests.cmake)
include(TES/Tests.cmake)
include(TwoPhaseFlowWithPP/Tests.cmake)
include(TwoPhaseFlowWithPrho/Tests.cmake)
