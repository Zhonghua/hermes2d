find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})
find_package(UMFPACK REQUIRED)
if(NOT UMFPACK_NO_BLAS)
	enable_language(Fortran)
	find_package(BLAS REQUIRED)
endif(NOT UMFPACK_NO_BLAS)

# examples
add_subdirectory(01-mesh)
add_subdirectory(02-space)
add_subdirectory(03-poisson)
add_subdirectory(04-bc-dirichlet)
add_subdirectory(05-bc-neumann)
add_subdirectory(06-bc-newton)
add_subdirectory(07-general)
add_subdirectory(08-system)
add_subdirectory(09-timedep)
add_subdirectory(10-adapt)
add_subdirectory(11-adapt-system)
add_subdirectory(12-adapt-general)
add_subdirectory(13-newton-elliptic-1)
add_subdirectory(14-newton-elliptic-2)
add_subdirectory(15-newton-elliptic-adapt)
add_subdirectory(16-newton-timedep-heat)
add_subdirectory(17-newton-timedep-flame)
add_subdirectory(18-newton-timedep-ns)
add_subdirectory(19-newton-timedep-gp)
add_subdirectory(20-newton-timedep-heat-adapt)
add_subdirectory(23-newton-timedep-gp-adapt)
