# Room-b118 Math Library
# signal - Sinal module

cmake_minimum_required(VERSION 3.23)
project(Room-b118-signal VERSION 1.0.0 LANGUAGES CXX)

add_library(b118-signal INTERFACE)

set_target_properties(b118-signal
    PROPERTIES
        # LINKER_LANGUAGE CXX
        CMAKE_CXX_STANDARD 17
        CMAKE_CXX_STANDARD_REQUIRED True
)

add_dependencies(b118-signal b118-core)

target_sources(b118-signal
    PUBLIC
        FILE_SET b118_signal_public_headers
            TYPE HEADERS
            BASE_DIRS include/
            FILES
            include/b118/signal/convolution.hpp
)

target_include_directories(b118-signal INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>  # <prefix>/include
)

target_link_libraries(b118-linalg
    INTERFACE
        b118-core
        ${LIB_OPENBLAS}
        ${LIB_FFTW3}
)

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-signal
    EXPORT "b118-signal-targets"
    FILE_SET b118_signal_public_headers
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"  # include
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
)

# add_subdirectory(tests)