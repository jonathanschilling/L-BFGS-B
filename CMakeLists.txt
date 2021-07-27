cmake_minimum_required(VERSION 3.10)

# set the project name
project(lbfgsb Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# shared library for L-BFGS-B
set(liblbfgsb_src "")
add_subdirectory (src)
add_library(lbfgsb SHARED ${liblbfgsb_src})

if(BLAS_FOUND)
    target_link_libraries(lbfgsb PRIVATE ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

if(LAPACK_FOUND)
    target_link_libraries(lbfgsb PRIVATE ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

# example executables
add_executable(driver1_f77 drivers/driver1.f)
add_executable(driver1_f90 drivers/driver1.f90)

add_executable(driver2_f77 drivers/driver2.f)
add_executable(driver2_f90 drivers/driver2.f90)

add_executable(driver3_f77 drivers/driver3.f)
add_executable(driver3_f90 drivers/driver3.f90)

target_link_libraries(driver1_f77 PRIVATE lbfgsb)
target_link_libraries(driver1_f90 PRIVATE lbfgsb)

target_link_libraries(driver2_f77 PRIVATE lbfgsb)
target_link_libraries(driver2_f90 PRIVATE lbfgsb)

target_link_libraries(driver3_f77 PRIVATE lbfgsb)
target_link_libraries(driver3_f90 PRIVATE lbfgsb)
