cmake_minimum_required(VERSION 3.28.3)

project(pca VERSION 0.1
        DESCRIPTION "First compile with CMake"
        LANGUAGES Fortran)

enable_language(Fortran)

add_library(util_mod ../../src/spline_functions/utils.f90)

# file(GLOB_RECURSE sources src/spline_functions/utils.f90)

add_executable(read_input_file_2 read_input_file.f90    )

target_link_libraries(read_input_file_2 PRIVATE util_mod)

install(TARGETS read_input_file_2
        RUNTIME DESTINATION "bin")
# INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/src DESTINATION ./bin)