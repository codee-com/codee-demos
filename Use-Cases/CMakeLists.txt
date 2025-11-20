cmake_minimum_required(VERSION 3.20)
project(Use-Cases LANGUAGES C Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -c -cpp -Wall -Wuninitialized -Wmaybe-uninitialized")
set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       -c -Wall -Wuninitialized -Wmaybe-uninitialized")

file(GLOB F90_SRCS "*.f90")
file(GLOB C_SRCS   "*.c")

foreach(src ${F90_SRCS})
    get_filename_component(prog ${src} NAME_WE)
    add_executable(${prog} ${src})
endforeach()

foreach(src ${C_SRCS})
    get_filename_component(prog ${src} NAME_WE)
    add_executable(${prog} ${src})
endforeach()
