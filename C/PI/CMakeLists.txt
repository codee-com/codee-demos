cmake_minimum_required(VERSION 3.10)
project(PI C)

find_package(OpenMP REQUIRED)

# Set parallel flag for Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -parallel")
# Set parallel flag for cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qpar")
endif()

add_executable(pi pi.c)
target_link_libraries(pi PRIVATE m OpenMP::OpenMP_C)

add_custom_target(run
    COMMAND pi 1000000000
    DEPENDS pi
)
