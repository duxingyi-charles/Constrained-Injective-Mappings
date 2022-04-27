cmake_minimum_required(VERSION 3.15.0)
project(TLC VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")


include(FetchContent)
include(cmake/Eigen3.cmake)
include(cmake/cli11.cmake)
include(cmake/nlohmann-json.cmake)
include(cmake/nlopt.cmake)



aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

add_library(arcOverlap STATIC ${SRC_LIST})
target_include_directories(arcOverlap PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(arcOverlap PUBLIC Eigen3::Eigen)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(arcOverlap PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(findInjective apps/TLC.cpp)
target_link_libraries(findInjective PUBLIC arcOverlap nlopt::nlopt)
set_target_properties(findInjective PROPERTIES LINKER_LANGUAGE CXX)



add_executable(arc_arc_intersection_test apps/arc_arc_intersection_test.cpp)
target_link_libraries(arc_arc_intersection_test PUBLIC arcOverlap)

add_executable(subdivide_polyArc_test apps/subdivide_polyArc_test.cpp)
target_link_libraries(subdivide_polyArc_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(cell_decompose_test apps/cell_decompose_test.cpp)
target_link_libraries(cell_decompose_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arrangement_test apps/arrangement_test.cpp)
target_link_libraries(arrangement_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arc_occupancy_test apps/arc_occupancy_test.cpp)
target_link_libraries(arc_occupancy_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arc_occupancy_gradient_test apps/arc_occupancy_gradient_test.cpp)
target_link_libraries(arc_occupancy_gradient_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arcOverlap_test apps/arcOverlap_test.cpp)
target_link_libraries(arcOverlap_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arcOverlap_Iso_test apps/arcOverlap_Iso_test.cpp)
target_link_libraries(arcOverlap_Iso_test PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arcOverlap_QN apps/ArcOverlap_QN.cpp)
target_link_libraries(arcOverlap_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

add_executable(arcOverlap_PN apps/ArcOverlap_PN.cpp)
target_link_libraries(arcOverlap_PN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(arcOverlap_Iso_QN apps/ArcOverlap_Iso_QN.cpp)
target_link_libraries(arcOverlap_Iso_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

add_executable(sTLC_Iso_QN "apps/sTLC_Iso_QN.cpp")
target_link_libraries(sTLC_Iso_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

add_executable(sTGC_QN "apps/sTGC_QN.cpp")
target_link_libraries(sTGC_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

add_executable(sTGC_PN "apps/sTGC_PN.cpp")
target_link_libraries(sTGC_PN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json)

add_executable(sTGC_3D_QN "apps/sTGC_3D_QN.cpp")
target_link_libraries(sTGC_3D_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

add_executable(SEA_QN "apps/SEA_QN.cpp")
target_link_libraries(SEA_QN PUBLIC arcOverlap CLI11::CLI11 nlohmann::json nlopt::nlopt)

# suitesparse

set(SUITESPARSE_INCLUDE_DIR /usr/local/include CACHE PATH "suitesparse include directory")
target_include_directories(arcOverlap_PN PUBLIC "${SUITESPARSE_INCLUDE_DIR}")
target_include_directories(sTGC_PN PUBLIC "${SUITESPARSE_INCLUDE_DIR}")

if (WIN32)
    # ------------------------------------------------------------------
    # Detect SuiteSparse libraries:
    # If not found automatically, set SuiteSparse_DIR in CMake to the
    # directory where SuiteSparse-config.cmake was installed.
    # ------------------------------------------------------------------
    find_package(SuiteSparse CONFIG REQUIRED)
    target_link_libraries(arcOverlap_PN PRIVATE SuiteSparse::cholmod)
    target_link_libraries(sTGC_PN PRIVATE SuiteSparse::cholmod)
else()
    set(SUITESPARSE_LIB_DIR /usr/local/lib CACHE PATH "suitesparse library directory")
    target_link_libraries(arcOverlap_PN PUBLIC -L"${SUITESPARSE_LIB_DIR}" cholmod)
    target_link_libraries(sTGC_PN PUBLIC -L"${SUITESPARSE_LIB_DIR}" cholmod)
endif()


# Math library
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(findInjective PUBLIC ${MATH_LIBRARY})
    target_link_libraries(arc_arc_intersection_test PUBLIC ${MATH_LIBRARY})
    target_link_libraries(subdivide_polyArc_test PUBLIC ${MATH_LIBRARY})
endif()