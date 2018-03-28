APPEND_SOURCE_FILES(SOURCES)

add_library(PipeNetwork ${SOURCES})
target_link_libraries(PipeNetwork PUBLIC ProcessLib)

# include(Tests.cmake)