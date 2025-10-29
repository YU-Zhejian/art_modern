# TODO: Add compile options --rtlib=compiler-rt 
# TODO: Build clang using -DLIBCXX_USE_COMPILER_RT=YES -DLIBCXXABI_USE_COMPILER_RT=YES -DLIBCXXABI_USE_LLVM_UNWINDER=YES 

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_LINKER ld.lld)
set(CMAKE_AR llvm-ar)
set(CMAKE_NM llvm-nm)
set(CMAKE_RANLIB llvm-ranlib)
add_compile_options("-stdlib=libc++")
add_link_options("-stdlib=libc++")
add_link_options("-fuse-ld=lld")
