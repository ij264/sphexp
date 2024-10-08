cmake_minimum_required(VERSION 3.5)

project(SphExp)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

find_package(LAPACK)

add_library(SphExpModules src/module_function.f90
			              src/module_util.f90
			              src/nrtype.f90
			              src/nrutil.f90
						  src/module_sphexp.f90)






add_executable(sphexp_inv src/sphexp_inv.f90)
target_link_libraries(sphexp_inv  SphExpModules lapack)





