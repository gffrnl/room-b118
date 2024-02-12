# Room b118 math library

cmake_minimum_required(VERSION 3.12)
project(
    b118
    VERSION 1.0.0
    LANGUAGES CXX
)

#enable_language(C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_BUILD_TYPE DEBUG)

add_library(
    b118
    include/b118/almost_equal.hpp
    include/b118/equal_within_ulps.hpp
    include/b118/linspace.hpp
    src/file.cpp # CMake can not determine linker language for target without a .cpp file

    # Vendors
    vendor/cppreference/include/b118/vendor/cppreference/equal_within_ulps.hpp
)
target_include_directories(
    b118
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/vendor/cppreference/include"
)
target_compile_definitions(
    b118
    PUBLIC B118_VERSION="${CMAKE_PROJECT_VERSION}"
    PUBLIC B118_VERSION_MAJOR="${CMAKE_PROJECT_VERSION_MAJOR}"
    PUBLIC B118_VERSION_MINOR="${CMAKE_PROJECT_VERSION_MINOR}"
    PUBLIC B118_VERSION_PATCH="${CMAKE_PROJECT_VERSION_PATCH}"
    PUBLIC B118_VERSION_TWEAK="${CMAKE_PROJECT_VERSION_TWEAK}"
)

add_subdirectory(examples)