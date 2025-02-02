project(slim_libfmt LANGUAGES CXX)
include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories("${CMAKE_CURRENT_LIST_DIR}/include")

set(CEU_CM_SILENT ON)
set(CEU_CM_CC_FLAGS "")
include("${LIBCEU_CMAKE_COLLECTION_DIR}/get_compiler_flags.cmake")
set(SLIM_LIBFMT_CC_FLAGS ${CEU_CM_CC_FLAGS})
unset(CEU_CM_SILENT)

file(GLOB SLIM_LIBFMT_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*/*.cc)
if(BUILD_SHARED_LIBS)
    add_library(slim_libfmt SHARED ${SLIM_LIBFMT_SOURCES})
else()
    add_library(slim_libfmt STATIC ${SLIM_LIBFMT_SOURCES})
endif()

message(STATUS "slim_libfmt compiler flags: ${SLIM_LIBFMT_CC_FLAGS}")
target_compile_options(slim_libfmt BEFORE PRIVATE "${SLIM_LIBFMT_CC_FLAGS}")
