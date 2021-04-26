cmake_minimum_required(VERSION 3.15.0)
project(TLC VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")


include(FetchContent)
include(cmake/Eigen3.cmake)
include(cmake/nlopt.cmake)



aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

add_library(arcOverlap STATIC ${SRC_LIST})
target_include_directories(arcOverlap PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(arcOverlap PUBLIC Eigen3::Eigen)

add_executable(findInjective apps/TLC.cpp)
target_link_libraries(findInjective PUBLIC arcOverlap nlopt::nlopt)
set_target_properties(findInjective PROPERTIES LINKER_LANGUAGE CXX)



add_executable(arc_arc_intersection_test apps/arc_arc_intersection_test.cpp)
target_link_libraries(arc_arc_intersection_test PUBLIC arcOverlap)


find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(findInjective PUBLIC ${MATH_LIBRARY})
    target_link_libraries(arc_arc_intersection_test PUBLIC ${MATH_LIBRARY})
endif()