if(WITH_PROTOBUF)
    find_package(Protobuf REQUIRED)

    include_directories(${Protobuf_INCLUDE_DIRS})

    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} protobuf::libprotobuf)
endif ()