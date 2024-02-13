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




set(B188_INCLUDE_DIR include)

add_subdirectory(${B188_INCLUDE_DIR})


message(STATUS "B118_MODULE_BASIC_HEADERS    : ${B118_MODULE_BASIC_HEADERS}")
message(STATUS "B118_MODULE_LINALG_HEADERS   : ${B118_MODULE_LINALG_HEADERS}")
message(STATUS "B118_MODULE_FRLAP_HEADERS    : ${B118_MODULE_FRLAP_HEADERS}")


# ##############################################################################
# Basic module
# ##############################################################################

add_library(
    b118-basic
    ${B118_MODULE_BASIC_HEADERS}
    src/file.cpp # CMake can not determine linker language for target without a .cpp file
    # Vendors
    vendor/cppreference/include/b118/vendor/cppreference/equal_within_ulps.hpp
)
target_include_directories(
    b118-basic
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/vendor/cppreference/include"
)

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)


foreach(header ${B118_MODULE_BASIC_HEADERS})
  file(RELATIVE_PATH header_file_path
       "${CMAKE_CURRENT_SOURCE_DIR}/include"
       "${CMAKE_CURRENT_SOURCE_DIR}/${header}")
  get_filename_component(header_directory_path "${header_file_path}" DIRECTORY)
  message(STATUS "${header_file_path} on directory ${header_directory_path}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${header}"
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${header_directory_path}")
endforeach()
#vendor headres
file(RELATIVE_PATH header_file_path
       "${CMAKE_CURRENT_SOURCE_DIR}/vendor/cppreference/include/"
       "${CMAKE_CURRENT_SOURCE_DIR}/vendor/cppreference/include/b118/vendor/cppreference/equal_within_ulps.hpp")
get_filename_component(header_directory_path "${header_file_path}" DIRECTORY)
message(STATUS "${header_file_path} on directory ${header_directory_path}")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/vendor/cppreference/include/b118/vendor/cppreference/equal_within_ulps.hpp"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${header_directory_path}")

# Debug Suffix
#set_target_properties(${library_name} PROPERTIES DEBUG_POSTFIX "d")

# Destinations

# if you haven't included it already
# definitions of CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR and others
#include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-basic
  EXPORT "b118-basicTargets"
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)



add_library(
    b118
    # Basic module
    ${B118_MODULE_BASIC_HEADERS}
    # Linear algebra module b118::linalg
    ${B118_MODULE_LINALG_HEADERS}
    # Fractional Laplacian module: b118::frlap
    ${B118_MODULE_FRLAP_HEADERS}

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