cmake_minimum_required(VERSION 2.8)

project(microbioreactor_parallel)
add_subdirectory (./external/Biser)
set (PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

find_package(MPI REQUIRED)

find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "-std=c++11 -O2 -Wall -pedantic -Wsign-compare ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

execute_process(
	COMMAND  git submodule update --init --recursive
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


set (PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set (SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set (INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

include_directories(SYSTEM ${MPI_INCLUDE_PATH} "${SRC_DIR}" "${INC_DIR}" "./external/Biser/include")

#file (GLOB_RECURSE SRC_FILES "${SRC_DIR}/*.cxx")

add_library(algs ${SRC_DIR}/algos.cxx)
target_link_libraries(algs ${MPI_C_LIBRARIES} calculator)

add_executable("sample3" ${SRC_DIR}/tmp3.cxx)
target_link_libraries("sample3" "${MPI_C_LIBRARIES}" "${OpenMP_LIBRARIES}" algs)

## Steady state model
add_executable("steady" ${SRC_DIR}/steady.cxx)
target_link_libraries("steady" "${MPI_C_LIBRARIES}" algs)

## Time changing model, reference, no parallel
add_executable("reference_t1" ${SRC_DIR}/reference_t1.cxx)
target_link_libraries("reference_t1" algs calculator)

## Parallel of time changing model, OpenMP
add_executable("parallel_omp" ${SRC_DIR}/parallel_omp.cxx)
target_link_libraries("parallel_omp" "${OpenMP_LIBRARIES}" algs calculator)

## Parallel of time changing model, MPI
add_executable("parallel_mpi" ${SRC_DIR}/parallel_mpi.cxx)
target_link_libraries("parallel_mpi" "${MPI_C_LIBRARIES}" algs calculator)
