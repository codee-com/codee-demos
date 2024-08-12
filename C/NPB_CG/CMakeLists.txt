cmake_minimum_required(VERSION 3.10)
project(NPB C)

set(BENCHMARK cg)
set(BENCHMARKU CG)
set(BENCHMARK_CLASS B)

#---------------------------------------------------------------------------
# The variable RAND controls which random number generator 
# is used. It is described in detail in README.install. 
# Use "randi8" unless there is a reason to use another one. 
# Other allowed values are "randi8_safe", "randdp" and "randdpvec"
#---------------------------------------------------------------------------
# set(RAND randi8)
# The following is highly reliable but may be slow:
set(RAND randdp)

#---------------------------------------------------------------------------
# The variable WTIME is the name of the wtime source code module in the
# common directory.  
# For most machines,       use wtime.c
# For SGI power challenge: use wtime_sgi64.c
#---------------------------------------------------------------------------
set(WTIME  wtime.c)

add_subdirectory(sys)
add_subdirectory(common)
add_subdirectory(CG)

add_custom_target(run
    COMMAND ${CMAKE_BINARY_DIR}/bin/${BENCHMARK}.${BENCHMARK_CLASS}.x
    DEPENDS ${BENCHMARK}.${BENCHMARK_CLASS}.x 
)
