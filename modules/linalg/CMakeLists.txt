# Room-b118 Math Library
# linalg - Numerical linear algebra module

cmake_minimum_required(VERSION 3.23)
project(Room-b118-linalg VERSION 1.0.0 LANGUAGES CXX)

add_library(b118-linalg INTERFACE)

set_target_properties(b118-linalg
    PROPERTIES
        # LINKER_LANGUAGE CXX
        CMAKE_CXX_STANDARD 17
        CMAKE_CXX_STANDARD_REQUIRED True
)

add_dependencies(b118-linalg b118-core)

target_sources(b118-linalg
    PUBLIC
        FILE_SET b118_linalg_public_headers
            TYPE HEADERS
            BASE_DIRS include/
            FILES
            include/b118/linalg/matrix_kind.hpp
            include/b118/linalg/matrix/matrix.hpp
            include/b118/linalg/matrix/column.hpp
            include/b118/linalg/matrix/symmetric.hpp
            include/b118/linalg/matrix/symmetric_toeplitz.hpp
            include/b118/linalg/matrix.hpp
            include/b118/linalg/vector.hpp
            include/b118/linalg/toeplitz/fast_symm_toeplitz_prod.hpp
            include/b118/linalg/toeplitz/solvers/durbin.hpp
            include/b118/linalg/toeplitz/solvers/levinson.hpp
            include/b118/linalg/toeplitz/solvers/trench.hpp
)

target_include_directories(b118-linalg INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>  # <prefix>/include
)

target_link_libraries(b118-linalg
    INTERFACE
        b118-core
        ${LIB_OPENBLAS}
        ${LIB_FFTW3}
)


# Debug Postffix
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(b118-linalg PROPERTIES DEBUG_POSTFIX "d")
endif()

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-linalg
    EXPORT "b118-linalg-targets"
    FILE_SET b118_linalg_public_headers
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"  # include
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
)

add_subdirectory(tests)