unset(USE_BTREE_MAP_QUALGEN)
unset(USE_STL_QUALGEN)

if(DEFINED USE_BTREE_MAP)
    message(WARNING "-DUSE_BTREE_MAP is deprecated. Ignored.")
endif()

if(NOT DEFINED USE_QUAL_GEN)
    set(USE_QUAL_GEN WALKER)
endif()

if("${USE_QUAL_GEN}" STREQUAL "WALKER")
    set(USE_WALKER_QUALGEN ON)
elseif("${USE_QUAL_GEN}" STREQUAL "STL")
    set(USE_STL_QUALGEN ON)
else()
    message(FATAL_ERROR "Unknown USE_QUAL_GEN: ${USE_QUAL_GEN}. Should be one of WALKER, and STL")
endif()
