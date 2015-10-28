if(OGS_USE_PETSC)

add_executable(allgather_benchmark
    allgather_benchmark.cpp
)

target_link_libraries(allgather_benchmark
    ApplicationsLib
    BaseLib
    MathLib
)

target_link_libraries(allgather_benchmark ${PETSC_LIBRARIES})

endif()
