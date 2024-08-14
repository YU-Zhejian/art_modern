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
    set(CEU_CM_PRINTED_${NAME}
        TRUE
        CACHE INTERNAL "")
    message(
        STATUS
            "CEU_CMAKE_CHECK: Checking ${NAME} ... shared: compile=${CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_SHARED}, run=${CEU_CM_HAVE_WORKING_${VARNAME}_RUN_SHARED}; static: compile=${CEU_CM_HAVE_WORKING_${VARNAME}_COMPILE_STATIC}, run=${CEU_CM_HAVE_WORKING_${VARNAME}_RUN_STATIC}"
    )
endmacro()
