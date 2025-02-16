# TODO: Refactor this sh*t.
include("${CMAKE_CURRENT_LIST_DIR}/../enable_debug.cmake")

if(NOT DEFINED CEU_CM_TRY_COMPILE_INCLUDES)
    set(CEU_CM_TRY_COMPILE_INCLUDES
        ""
        CACHE INTERNAL "Extra path added to try_compile")
endif()

function(ceu_cm_enhanced_try_run)
    # Parse args
    set(options STATIC)
    set(oneValueArgs SRC_PATH VARNAME)
    set(multiValueArgs LINK_LIBRARIES LINK_FLAGS COMPILE_DEFS DEPENDS)
    cmake_parse_arguments(CEU_CM_ENHANCED_TRY_RUN "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(CEU_CM_ENHANCED_TRY_RUN_STATIC)
        set(TARGET_POSTFIX "STATIC")
    else()
        set(TARGET_POSTFIX "SHARED")
    endif()

    # Unconditional skip
    if(DEFINED SKIP_ALL_TRY_RUN)
        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
            ON
            CACHE INTERNAL "SKIP_ALL_TRY_RUN was set, skip all tests.")
        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
            0
            CACHE INTERNAL "SKIP_ALL_TRY_RUN was set, skip all tests.")
        return()
    endif()

    # Start processing
    if(NOT DEFINED CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
       OR NOT CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX} EQUAL 0)

        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
            OFF
            CACHE INTERNAL "Default to fail.")
        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
            127
            CACHE INTERNAL "Default to fail.")

        # Finding whether all dependencies are defined.
        if(DEFINED CEU_CM_ENHANCED_TRY_RUN_DEPENDS)
            foreach(DEPENDENT_VARNAME ${CEU_CM_ENHANCED_TRY_RUN_DEPENDS})

                # Assess compiling result first.
                if(DEFINED CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_COMPILE_${TARGET_POSTFIX})
                    if(NOT ${CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_COMPILE_${TARGET_POSTFIX}})
                        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
                            OFF
                            CACHE
                                INTERNAL
                                "Dependency CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_COMPILE_${TARGET_POSTFIX} failed.")
                    endif()
                else()
                    message(
                        FATAL_ERROR
                            "Dependency CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_COMPILE_${TARGET_POSTFIX} is not undefined"
                    )
                endif()

                # Assess running result.
                if(DEFINED CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_RUN_${TARGET_POSTFIX})
                    if(NOT ${CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_RUN_${TARGET_POSTFIX}} EQUAL 0)
                        set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
                            127
                            CACHE INTERNAL
                                  "Dependency CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_RUN_${TARGET_POSTFIX} failed.")
                        return()
                    endif()
                else()
                    set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
                        OFF
                        CACHE INTERNAL
                              "Dependency CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_RUN_${TARGET_POSTFIX} undefined.")
                    message(
                        FATAL_ERROR
                            "Dependency CEU_CM_HAVE_WORKING_${DEPENDENT_VARNAME}_RUN_${TARGET_POSTFIX} is not undefined"
                    )
                endif()
            endforeach()
        endif()

        if(DEFINED CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARIES)
            foreach(CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARY ${CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARIES})
                if(CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARY EQUAL "") # If the library is empty
                    continue()
                elseif(NOT CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARY) # If the library is not found
                    set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
                        OFF
                        CACHE INTERNAL "Dependency ${CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARY} not found.")
                    set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
                        127
                        CACHE INTERNAL "Dependency ${CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARY} not found.")
                    return()
                endif()
            endforeach()
        else()
            set(CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARIES "")
        endif()

        set(THIS_COMPILE_DEFS "")

        if(DEFINED CEU_CM_ENHANCED_TRY_RUN_COMPILE_DEFS)
            foreach(FLAG ${CEU_CM_ENHANCED_TRY_RUN_COMPILE_DEFS})
                set(THIS_COMPILE_DEFS ${THIS_COMPILE_DEFS} -D${FLAG})
            endforeach()
        endif()

        if(DEFINED CEU_CM_ADDITIONAL_COMPILER_FLAGS)
            foreach(FLAG ${CEU_CM_ADDITIONAL_COMPILER_FLAGS})
                set(THIS_COMPILE_DEFS ${THIS_COMPILE_DEFS} ${FLAG})
            endforeach()
        endif()

        if(CEU_CM_ENHANCED_TRY_RUN_STATIC AND NOT BORLAND) # Borland compiler does not support such functions
            set(THIS_LINK_FLAGS -static -static-libgcc -static-libstdc++ -static-libgfortran)
        else()
            set(THIS_LINK_FLAGS "")
        endif()

        if(DEFINED CEU_CM_ENHANCED_TRY_RUN_LINK_FLAGS)
            foreach(FLAG ${CEU_CM_ENHANCED_TRY_RUN_LINK_FLAGS})
                set(THIS_LINK_FLAGS ${THIS_LINK_FLAGS} ${FLAG})
            endforeach()
        endif()

        try_run(
            CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
            CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
            "${CMAKE_BINARY_DIR}/CEU_TRC" "${CEU_CM_ENHANCED_TRY_RUN_SRC_PATH}"
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CEU_CM_TRY_COMPILE_INCLUDES} -DLINK_FLAGS=${THIS_LINK_FLAGS}"
                        LINK_LIBRARIES ${CEU_CM_ENHANCED_TRY_RUN_LINK_LIBRARIES}
            COMPILE_DEFINITIONS ${THIS_COMPILE_DEFS}
            COMPILE_OUTPUT_VARIABLE CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}_VAR
            RUN_OUTPUT_VARIABLE CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}_VAR)

        if(NOT DEFINED CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}_VAR)
            set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}_VAR "")
        endif()
        if(NOT DEFINED CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}
           OR NOT CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX})
            file(
                WRITE
                "${CMAKE_BINARY_DIR}/compile_logs/test_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_${TARGET_POSTFIX}_compile.log"
                ${CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}_VAR})
        endif()

        unset(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_COMPILE_${TARGET_POSTFIX}_VAR)
        unset(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}_VAR)

        if(NOT DEFINED CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX})
            set(CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}
                127
                CACHE INTERNAL "Compilation failed.")
        endif()
        if(NOT CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX} EQUAL 0)
            file(WRITE
                 "${CMAKE_BINARY_DIR}/compile_logs/test_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_${TARGET_POSTFIX}_run.log"
                 ${CEU_CM_HAVE_WORKING_${CEU_CM_ENHANCED_TRY_RUN_VARNAME}_RUN_${TARGET_POSTFIX}_VAR})
        endif()
    else()
        message(DEBUG "MISC_ENHANCED_TRY_COMPILE: Cached static ${CEU_CM_ENHANCED_TRY_RUN_VARNAME}")
    endif()
endfunction()
