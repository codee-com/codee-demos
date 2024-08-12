cmake_minimum_required(VERSION 3.10)
project(CANNY C)

find_package(OpenMP REQUIRED)

# Set parallel flag for Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -parallel")
# Set parallel flag for cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qpar")
endif()

set(ZIP_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../15360_8640.zip")

include_directories(include)

add_executable(canny canny.c)
target_link_libraries(canny PRIVATE m OpenMP::OpenMP_C)

add_custom_target(testvecs
    COMMAND ${CMAKE_COMMAND} -E tar xf ${ZIP_FILE}
    DEPENDS ${ZIP_FILE}
)

add_custom_target(run
    COMMAND canny testvecs/input/15360_8640.pgm 0.5 0.7 0.9
    DEPENDS testvecs canny
)
