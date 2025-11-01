cmake_minimum_required(VERSION 3.28.3)

project(pca
        VERSION 0.1
        DESCRIPTION "First compile with CMake"
        LANGUAGES Fortran)

enable_language(Fortran)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=f2008")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} -fcheck=all -Wall -g -fbounds-check -Wno-uninitialized")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/src)

# Source code
add_subdirectory(src)

# Examples code
add_subdirectory(examples)