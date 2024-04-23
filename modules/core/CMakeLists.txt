# Room-b118 Math Library
# core - core module

cmake_minimum_required(VERSION 3.23)
project(Room-b118-core VERSION 1.0.0 LANGUAGES CXX)

add_library(b118-core INTERFACE)

set_target_properties(b118-core
    PROPERTIES
        # LINKER_LANGUAGE CXX
        CMAKE_CXX_STANDARD 17
        CMAKE_CXX_STANDARD_REQUIRED True
)

target_sources(b118-core
    PUBLIC
        FILE_SET b118_core_public_headers
            TYPE HEADERS
            BASE_DIRS include/
            FILES
                include/b118/utility.hpp
                include/b118/type_traits.hpp
                include/b118/numbers.hpp
                include/b118/almost_equal.hpp
                include/b118/equal_within_ulps.hpp
                include/b118/closest_sorted.hpp
                include/b118/linspace.hpp
                include/b118/grid.hpp
                include/b118/grid_function.hpp
                include/b118/fftw.hpp
                # vendors
                include/b118/vendor/cppreference/equal_within_ulps.hpp
)

target_include_directories(b118-core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>  # <prefix>/include
)

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-core
    EXPORT "b118-core-targets"
    FILE_SET b118_core_public_headers
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"  # include/b118
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
)
