cmake_minimum_required(VERSION 3.16)
project(tc2038_hull_lab CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)   # timings are meaningless without optimizations
endif()
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra)   # la rúbrica pide compilar sin warnings
endif()

# --- Partes 1, 2 y 5: no necesitan CGAL ---
add_library(hull hull.cpp)
add_executable(main main.cpp)            # banco de pruebas del equipo
target_link_libraries(main hull)
target_compile_definitions(main PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
add_executable(hull_main main_hull.cpp)
target_link_libraries(hull_main hull)
add_executable(app app.cpp)
target_link_libraries(app hull)

# --- Partes 0, 3, 4, 6 y bonus: requieren CGAL ---
find_package(CGAL QUIET)
if(NOT CGAL_FOUND)
  message(WARNING "CGAL no encontrado: solo se compilan main, hull_main y app")
  return()
endif()

add_executable(cgal_hello cgal_hello.cpp)
target_link_libraries(cgal_hello CGAL::CGAL)
add_executable(stress stress.cpp)
target_link_libraries(stress hull CGAL::CGAL)
add_executable(bench bench.cpp)
target_link_libraries(bench hull CGAL::CGAL)
add_executable(robust robust.cpp)
target_link_libraries(robust CGAL::CGAL)
add_executable(hull3d hull3d.cpp)
target_link_libraries(hull3d CGAL::CGAL)
