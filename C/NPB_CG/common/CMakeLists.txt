cmake_minimum_required(VERSION 3.10)
project(common C)

find_package(OpenMP REQUIRED)

# Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -parallel -O3 -g -Wall -mcmodel=medium")
# Intel icx
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -g -Wall -mcmodel=medium")
# GNU gcc
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -g -Wall -mcmodel=medium")
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
# LLVM clang-cl
  if (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} /O2 /W3")
# LLVM clang
    elseif (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -g -Wall -mcmodel=medium")
  endif()
# cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} /O2 /Qpar /W3")
# Other compilers
else()
    message("** Warning, compiler '${CMAKE_C_COMPILER_ID}' is not covered.")
    message("** It is possible that the flags are not the correct ones")
endif()

set (BINDIR	${CMAKE_CURRENT_SOURCE_DIR}/../bin)
set (BENCHMARK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../${BENCHMARKU})
set(COMMON ${CMAKE_CURRENT_SOURCE_DIR})


# Class "U" is used internally by the setparams program to mean
# "unknown". This means that if you don't specify CLASS=
# on the command line, you'll get an error. It would be nice
# to be able to avoid this, but we'd have to get information
# from the setparams back to the make program, which isn't easy.
set(BENCHMARK_CLASS "U")

# set (PROGRAM  $(BINDIR)/$(BENCHMARK).$(BENCHMARK_CLASS).x)

add_library(print_results OBJECT print_results.c)
add_library(c_print_results OBJECT c_print_results.c)
add_library(randdp OBJECT randdp.c)
add_library(c_timers OBJECT c_timers.c)
add_library(wtime OBJECT wtime.c)
