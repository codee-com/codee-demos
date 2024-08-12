cmake_minimum_required(VERSION 3.10)
project(MATMUL Fortran)

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    message("OpenMP is enabled")
else()
    message("OpenMP is disabled")
endif()

add_executable(matmul
    matmul.f90
    time.f90
    main.f90
)

add_custom_target(run
    COMMAND matmul 1500
    DEPENDS matmul
)
