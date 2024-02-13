# Room b118 math library

cmake_minimum_required(VERSION 3.23)
project(
    b118
    VERSION 1.0.0
    LANGUAGES CXX
)

#enable_language(C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_BUILD_TYPE DEBUG)


set(B188_INCLUDE_DIR include)
set(B188_SRC_DIR src)


add_subdirectory(${B188_INCLUDE_DIR})
add_subdirectory(${B188_SRC_DIR})

#include(module-basic.cmake)
include(module-frlap.cmake)

add_library(
    b118
    # Basic module
    ${B118_MODULE_BASIC_PUBLIC_HEADERS}
    # Linear algebra module b118::linalg
    ${B118_MODULE_LINALG_HEADERS}
    # Fractional Laplacian module: b118::frlap
    ${B118_MODULE_FRLAP_HEADERS}

    src/file.cpp # CMake can not determine linker language for target without a .cpp file
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