#[=======================================================================[
ceu_cm_print_test_status -- Print pretty-formatted test status.

Synopsis: ceu_cm_print_test_status(NAME VARNAME)

Params:
    - `NAME`: Name of the test
    - `VARNAME`: See below.

Requires:
    - `CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_SHARED`
    - `CEU_CM_HAVE_WORKING_${VARNAME}_RUN_SHARED`
    - `CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_STATIC`
    - `CEU_CM_HAVE_WORKING_${VARNAME}_RUN_STATIC`
#]=======================================================================]
macro(ceu_cm_print_test_status NAME VARNAME)
    if(DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_SHARED
       AND DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_RUN_SHARED
       AND DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_STATIC
       AND DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_RUN_STATIC)
        message(
            STATUS
                "CEU_CM_CHECK: Checking ${NAME} ... shared: compile=${CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_SHARED}, run=${CEU_CM_HAVE_WORKING_${VARNAME}_RUN_SHARED}"
        )
    endif()
    if(DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_STATIC AND DEFINED CEU_CM_HAVE_WORKING_${VARNAME}_RUN_STATIC)
        message(
            STATUS
                "CEU_CM_CHECK: Checking ${NAME} ... static: compile=${CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_STATIC}, run=${CEU_CM_HAVE_WORKING_${VARNAME}_RUN_STATIC}"
        )
    endif()
endmacro()
