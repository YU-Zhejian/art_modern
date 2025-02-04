#[=======================================================================[
ceu_cm_set_static_target -- Mark a target as static.

Synopsis: ceu_cm_set_static_target(name)

Params:
    - `name`: the name of the target that would be sattis

Notice:
    This function works by adding following lines to link options:
    - `-static`
    - `-static-libgcc`
    - `-static-libstdc++`

Warnings:
    - For those who uses pthreads & openMP on GLibC platforms (i.e., common GNU/Linux EXCLUDING Alpine Linux), this option would lead to segfaults.
#]=======================================================================]
macro(ceu_cm_set_static_target name)
    set_target_properties("${name}" PROPERTIES LINK_SEARCH_START_STATIC 1)
    set_target_properties("${name}" PROPERTIES LINK_SEARCH_END_STATIC 1)
    set_target_properties("${name}" PROPERTIES INSTALL_RPATH "")
    if(NOT (BORLAND OR MSVC))
        set_target_properties("${name}" PROPERTIES LINK_FLAGS -static $<$<COMPILE_LANGUAGE:C,CXX>:-static-libgcc>
                                                                      $<$<COMPILE_LANGUAGE:CXX>:-static-libstdc++>)
    endif()
endmacro()
