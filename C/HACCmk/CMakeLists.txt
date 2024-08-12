cmake_minimum_required(VERSION 3.10)
project(HACCmk C)


find_package(OpenMP REQUIRED)

# Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -parallel")
    set(CMAKE_C_FLAGS_ "-O0")
    set(CMAKE_C_FLAGS_MAIN "-g -O0 ${OpenMP_C_FLAGS} -Wall")
    set(CMAKE_C_FLAGS_STEP10 "-g -O3 -xAVX2 -Wall")
# Intel icx
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(CMAKE_C_FLAGS_ "-O0")
    set(CMAKE_C_FLAGS_MAIN "-g -O0 ${OpenMP_C_FLAGS} -Wall")
    set(CMAKE_C_FLAGS_STEP10 "-g -O3 -xAVX2 -Wall")
# GNU gcc
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS_ "-O0")
    set(CMAKE_C_FLAGS_MAIN "-g -O0 ${OpenMP_C_FLAGS} -Wall")
    set(CMAKE_C_FLAGS_STEP10 "-g -O3 -mavx2 -Wall")
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
# LLVM clang-cl
  if (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "MSVC")
    set(CMAKE_C_FLAGS_ "/Od")
    set(CMAKE_C_FLAGS_MAIN "/Od ${OpenMP_C_FLAGS} /W3")
    set(CMAKE_C_FLAGS_STEP10 "/O2 /arch:AVX2 /W3")
# LLVM clang
    elseif (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "GNU")
        set(CMAKE_C_FLAGS_ "-O0")
        set(CMAKE_C_FLAGS_MAIN "-g -O0 ${OpenMP_C_FLAGS} -Wall")
        set(CMAKE_C_FLAGS_STEP10 "-g -O3 -mavx2 -Wall")
  endif()
# cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qpar")
    set(CMAKE_C_FLAGS_ "/Od")
    set(CMAKE_C_FLAGS_MAIN "/Od ${OpenMP_C_FLAGS} /W3")
    set(CMAKE_C_FLAGS_STEP10 "/O2 /arch:AVX2 /W3")
# Other compilers
else()
    message("** Warning, compiler '${CMAKE_C_COMPILER_ID}' is not covered.")
    message("** It is possible that the flags are not the correct ones")
    set(CMAKE_C_FLAGS_ "-O0")
    set(CMAKE_C_FLAGS_MAIN "-g -O0 ${OpenMP_C_FLAGS} -Wall")
    set(CMAKE_C_FLAGS_STEP10 "-g -O3 -mavx2 -Wall")
endif()

set_source_files_properties(main.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS_MAIN}")
set_source_files_properties(Step10_orig.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS_STEP10}")
set_source_files_properties(mysecond.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS_}")

add_executable(main main.c Step10_orig.c mysecond.c)
target_link_libraries(main -L$ENV{HOME}/lib -lm)
target_link_options(main PRIVATE ${CMAKE_C_FLAGS_} ${OpenMP_C_FLAGS})

add_custom_target(run
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/main
    DEPENDS main
)
