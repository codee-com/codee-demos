cmake_minimum_required(VERSION 3.10)
project(Himeno Fortran)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    message("OpenMP is enabled")
else()
    message("OpenMP is disabled")
endif()

add_executable(himeno
    himeno.f90
)

add_custom_target(run
    COMMAND himeno 1500
    DEPENDS himeno
)

