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


#set(NLopt_INCLUDE_DIR /usr/local/include CACHE PATH "NLopt include directory")
#set(NLopt_LIB_DIR /usr/local/lib CACHE PATH "NLopt library directory")
#
#target_include_directories(findInjective PUBLIC "${NLopt_INCLUDE_DIR}")
#
#target_link_libraries(findInjective PUBLIC -L"${NLopt_LIB_DIR}" nlopt)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(findInjective PUBLIC ${MATH_LIBRARY})
endif()

