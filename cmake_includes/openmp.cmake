find_package(OpenMP)
if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
    include_directories(${OpenMP_C_INCLUDE_DIRS})
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    set(ART_MODERN_LINK_LIBS ${ART_MODERN_LINK_LIBS} OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
    set(WITH_OPENMP ON)
endif()

# OpenMP_<lang>_VERSION OpenMP version implemented by the <lang> compiler.
