APPEND_SOURCE_FILES(SOURCES)

add_library(TH2M ${SOURCES})
target_link_libraries(TH2M PUBLIC ProcessLib)

include(Tests.cmake)
