project(slim_libceu LANGUAGES CXX)
include_directories(${CMAKE_CURRENT_LIST_DIR})
include_directories("${CMAKE_BINARY_DIR}/include")

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)

check_include_file("stdnoreturn.h" CEU_HAVE_INCLUDE_STDNORETURN_H)
check_include_file("cygwin/version.h" CEU_HAVE_INCLUDE_CYGWIN_VERSION_H)
check_include_file("sys/utsname.h" CEU_HAVE_INCLUDE_SYS_UTSNAME_H)
check_include_file("_mingw.h" CEU_HAVE_INCLUDE__MINGW_H)

configure_file("${CMAKE_CURRENT_LIST_DIR}/libceu_stddef.h.in" "${CMAKE_BINARY_DIR}/include/libceu_stddef.h")

file(GLOB CEU_LIBRARY_SOURCES ${CMAKE_CURRENT_LIST_DIR}/*/*.cc)
if(BUILD_SHARED_LIBS)
    add_library(slim_libceu SHARED ${CEU_LIBRARY_SOURCES})
else()
    add_library(slim_libceu STATIC ${CEU_LIBRARY_SOURCES})
endif()
install(TARGETS slim_libceu)

