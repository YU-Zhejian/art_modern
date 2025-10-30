# This is a CMake toolchain file for using the LLVM toolchain (Clang, LLD, etc.)
# Note that it still contain GNU components like libgcc/libatomic.
# For removing GNU components, consider using a full LLVM-based toolchain.
# However, you need to build almost everything from source for that.
# I figured it out using ~10 hours of work under the direction of <https://gist.github.com/PancakeTAS/949f52ecac1a0d202059a6e0c3c03b02>.
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_LINKER ld.lld)
set(CMAKE_AR llvm-ar)
set(CMAKE_NM llvm-nm)
set(CMAKE_RANLIB llvm-ranlib)
add_compile_options("-stdlib=libc++")
add_link_options("-stdlib=libc++")
add_link_options("-fuse-ld=lld")
