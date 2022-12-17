cmake_minimum_required(VERSION 3.23)
project(delaunay_triangulation_implementation)

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)

add_executable(
        delaunay_triangulation_implementation delaunay_triangulation.cpp
)

target_link_libraries(
        delaunay_triangulation_implementation
        CGAL::CGAL
)
