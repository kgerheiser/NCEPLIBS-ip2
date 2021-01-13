#list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

execute_process(COMMAND cmake -E create_symlink
  "${CMAKE_CURRENT_SOURCE_DIR}/data"
  "${CMAKE_CURRENT_BINARY_DIR}/data" # New name
  )

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-r8 -g -traceback -check all -warn all -assume byterecl")
  set(CMAKE_C_FLAGS "-std=c99")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fno-range-check -g -fbacktrace -fcheck=all -Wall -O0 -fimplicit-none -Wsurprising -Wextra")
endif()

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "^(GNU)$" AND ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch -fallow-invalid-boz")
endif()

# Contains common routines used for testing such as reading and interpolating data
add_library(test_library input_data_mod.f90 interp_mod.f90)
target_link_libraries(test_library PUBLIC ip2::ip2_d)

add_executable(test_gdswzd test_gdswzd.c)
add_executable(test_scalar test_scalar.f90)
add_executable(test_vector test_vector.f90)

target_link_libraries(test_gdswzd PRIVATE test_library)
target_link_libraries(test_scalar PRIVATE test_library)
target_link_libraries(test_vector PRIVATE test_library)

#add_test(test_gdswzd_c test_gdswzd)

# scalar tests
add_test(test_lambert_bilinear_scalar test_scalar 218 0)
add_test(test_gaussian_neighbor_scalar test_scalar 127 2)
add_test(test_latlon_bilinear_scalar test_scalar 3 0)
add_test(test_mercator_bicubic_scalar test_scalar 8 1)
#add_test(test_polar-stereo_neighbor-budget_scalar test_scalar 212 6)
#add_test(test_rotatedB_spectral_scalar test_scalar 205 4)
add_test(test_rotatedE_budget_scalar test_scalar 203 3)

# # vector tests
# add_test(test_lambert_biliner_vector test_vector 218 0)
# add_test(test_gaussian_neighbor_vector test_vector 127 2)
# add_test(test_latlon_bilinear_vector test_vector 3 0)
# add_test(test_mercator_bicubic_vector test_vector 8 1)
# add_test(test_polar-stereo_neighbor-budget_vector test_vector 212 6)
# add_test(test_rotatedB_spectral_vector test_vector 205 4)
# add_test(test_rotatedE_budget_vector test_vector 203 3)


add_library(test_library_grib1 input_data_mod_grib1.f90 interp_mod_grib1.f90)
target_link_libraries(test_library_grib1 PUBLIC ip2::ip2_d)

add_executable(test_gdswzd_grib1 test_gdswzd_grib1.c)
add_executable(test_scalar_grib1 test_scalar_grib1.f90)
#add_executable(test_vector_grib1 test_vector_grib1.f90)


target_link_libraries(test_scalar_grib1 PRIVATE test_library_grib1)
#target_link_libraries(test_vector_grib1 PRIVATE test_library_grib1)
target_link_libraries(test_gdswzd_grib1 ip2::ip2_d)

add_test(test_gdswzd_c_grib1 test_gdswzd_grib1)
add_test(test_lambert_bilinear_scalar_grib1 test_scalar_grib1 218 0)
add_test(test_gaussian_neighbor_scalar_grib1 test_scalar_grib1 127 2)
add_test(test_latlon_bilinear_scalar_grib1 test_scalar_grib1 3 0)
add_test(test_mercator_bicubic_scalar_grib1 test_scalar_grib1 8 1)
#add_test(test_polar-stereo_neighbor-budget_scalar test_scalar 212 6)
# add_test(test_rotatedB_spectral_scalar test_scalar 205 4)
add_test(test_rotatedE_budget_scalar_grib1 test_scalar_grib1 203 3)

#add_executable(test_scalar_grib1 test_scalar_grib1.f90)
#add_executable(test_vector_grib1 test_vector_grib1.f90)

#target_link_libraries(test_gdswzd_grib1 test_library_grib1)
# target_link_libraries(test_scalar_grib1 test_scalar_grib1)
# target_link_libraries(test_vector_grib1 test_vector_grib1)
