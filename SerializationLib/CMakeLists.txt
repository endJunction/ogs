
if(Protobuf_FOUND)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS common.proto)
    list(APPEND SOURCES ${PROTO_SRCS} ${PROTO_HDRS})
    protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS integration_point.proto)
    list(APPEND SOURCES ${PROTO_SRCS} ${PROTO_HDRS})
    protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS material.proto)
    list(APPEND SOURCES ${PROTO_SRCS} ${PROTO_HDRS})

add_library(SerializationLib ${SOURCES})
endif() # Protobuf_FOUND

if(Protobuf_FOUND)
    target_link_libraries(SerializationLib ${PROTOBUF_LIBRARIES})
endif() # Protobuf_FOUND
