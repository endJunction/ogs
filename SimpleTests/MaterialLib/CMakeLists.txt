add_executable(fluid_switch_vs_fctptr
    fluid_switch_vs_fctptr.cpp
    ${SOURCES}
    ${HEADERS}
)
add_dependencies(fluid_switch_vs_fctptr celero)
target_link_libraries(fluid_switch_vs_fctptr
    MaterialLib
    celero
)
