cmake_minimum_required(VERSION 3.10)
project(CG C)

find_package(OpenMP REQUIRED)

set(TARGET ${BENCHMARK}.${BENCHMARK_CLASS}.x)
set(COMMON ${CMAKE_CURRENT_SOURCE_DIR}/../common)
set(config ${CMAKE_CURRENT_SOURCE_DIR}/../config)
set(config_bin ${CMAKE_CURRENT_BINARY_DIR}/../config)
set(sys ${CMAKE_CURRENT_SOURCE_DIR}/../sys)
set(sys_bin ${CMAKE_CURRENT_BINARY_DIR}/../sys)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Intel icc
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -parallel ${OpenMP_C_FLAGS} -mcmodel=medium")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "-I${COMMON} -g -Wall -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
# Intel icx
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "-I${COMMON} -g -Wall -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
# GNU gcc
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "-I${COMMON} -g -Wall -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
# LLVM clang-cl
  if (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2 ${OpenMP_C_FLAGS}")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "/I${COMMON} /O2 ${OpenMP_C_FLAGS} /W3")
# LLVM clang
    elseif (${CMAKE_C_COMPILER_FRONTEND_VARIANT} STREQUAL "GNU")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
        set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "-I${COMMON} -g -Wall -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
  endif()
# cl.exe
elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2 ${OpenMP_C_FLAGS} /Qpar")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "/I${COMMON} /O2 ${OpenMP_C_FLAGS} /Qpar /W3")
# Other compilers
else()
    message("** Warning, compiler '${CMAKE_C_COMPILER_ID}' is not covered.")
    message("** It is possible that the flags are not the correct ones")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
    set_source_files_properties(${BENCHMARK}.c PROPERTIES COMPILE_FLAGS "-I${COMMON} -g -Wall -O3 ${OpenMP_C_FLAGS} -mcmodel=medium")
endif()

configure_file(${config}/make.def ${config_bin}/make.def COPYONLY)

# Generate npbparams.h based
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/npbparams.h
    COMMAND ${sys_bin}/setparams ${BENCHMARK} ${BENCHMARK_CLASS}
    DEPENDS setparams ${config}/make.def
)

add_library(${BENCHMARK} OBJECT ${BENCHMARK}.c npbparams.h globals.h)
target_include_directories(${BENCHMARK} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${BENCHMARK} PRIVATE ${sys})

add_executable(${TARGET}
    ${BENCHMARK}.c
    ${COMMON}/print_results.c
    ${COMMON}/${RAND}.c
    ${COMMON}/c_timers.c
    ${COMMON}/${WTIME}
)

target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(${TARGET} ${BENCHMARK})
target_link_libraries(${TARGET} PRIVATE m)
