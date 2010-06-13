project(convert_quads_to_triangles)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(convert_quads_to_triangles-1 "${BIN}" square.mesh 10000)
add_test(convert_quads_to_triangles-2 "${BIN}" square_tri.mesh 30000)
add_test(convert_quads_to_triangles-3 "${BIN}" domain.mesh 50000)
add_test(convert_quads_to_triangles-4 "${BIN}" bracket.mesh 70000)

