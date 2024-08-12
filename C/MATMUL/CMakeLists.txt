cmake_minimum_required(VERSION 3.10)
project(MATMUL C)

find_package(OpenMP REQUIRED)

# Set parallel flag for Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -parallel")
# Set parallel flag for cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qpar")
endif()

include_directories(include)

add_executable(matmul
    matrix.c
    clock.c
    main.c
)
target_link_libraries(matmul PRIVATE OpenMP::OpenMP_C)

add_custom_target(run
    COMMAND matmul 1500
    DEPENDS matmul
)
