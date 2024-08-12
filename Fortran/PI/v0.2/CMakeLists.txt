CMAKE_minimum_required( VERSION 3.5 )
enable_language( Fortran )
project( pi )

SET(CMAKE_Fortran_COMPILER gfortran)

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules )

file( GLOB pi_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 )
file( GLOB pi_MODULES ${CMAKE_Fortran_MODULE_DIRECTORY}/*.f90 )

include_directories (${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable( 
    pi
    ${pi_MODULES}
    ${pi_SOURCES} 
    )
