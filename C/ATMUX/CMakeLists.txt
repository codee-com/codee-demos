cmake_minimum_required(VERSION 3.10)
project(ATMUX C)

find_package(OpenMP REQUIRED)

# Set parallel flag for Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -parallel")
# Set parallel flag for cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qpar")
endif()

include_directories(lib)

add_executable(atmux
    lib/Matrix2D.c
    lib/Vector.c
    lib/CRSMatrix.c
    atmux.c
)
target_link_libraries(atmux PRIVATE OpenMP::OpenMP_C)
set_property(TARGET atmux PROPERTY C_STANDARD 99)

add_custom_target(run
    COMMAND atmux 10000
    DEPENDS atmux
)
