# Room-b118 Math Library
# frlap - Fractional Laplacian module

cmake_minimum_required(VERSION 3.23)
project(Room-b118-frlap VERSION 1.0.0 LANGUAGES CXX)

add_library(b118-frlap INTERFACE)

set_target_properties(b118-frlap
    PROPERTIES
    # LINKER_LANGUAGE CXX
    CMAKE_CXX_STANDARD 17
    CMAKE_CXX_STANDARD_REQUIRED True
    )

add_dependencies(b118-frlap b118-core)
add_dependencies(b118-frlap b118-signal)

target_sources(b118-frlap
    PUBLIC
        FILE_SET b118_frlap_public_headers
            TYPE HEADERS
            BASE_DIRS include/
            FILES
                include/b118/frlap.hpp
                include/b118/frlap/gdm.hpp
                # include/b118/frlap/gdm/strategy.hpp
                include/b118/frlap/gdm/trunc_uniform.hpp
                # gdm coefficients
                include/b118/frlap/gdm/coefficients/generator.hpp
                include/b118/frlap/gdm/coefficients/gorenflo_mainardi.hpp
                include/b118/frlap/gdm/coefficients/huang_oberman_1.hpp
                include/b118/frlap/gdm/coefficients/huang_oberman_2.hpp
                include/b118/frlap/gdm/coefficients/centered_3_point_periodized.hpp
                # gdm
                include/b118/frlap/gdm/truncated_uniform.hpp
                include/b118/frlap/gdm/far_field_estimator.hpp
                # # gdm strategies
                # include/b118/frlap/gdm/strategies/spectral.hpp
                # include/b118/frlap/gdm/strategies/spectral_qawo.hpp
                # include/b118/frlap/gdm/strategies/spectral_tanh_sinh.hpp
                # include/b118/frlap/gdm/strategies/gorenflo_mainardi.hpp
                # include/b118/frlap/gdm/strategies/huang_oberman_linear.hpp
                # include/b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp
                # include/b118/frlap/gdm/strategies/centered_3_point_periodized.hpp
)

# target_sources(b118-frlap
#     PRIVATE
#         src/b118-frlap.cpp
#         # gdm strategies
#         src/gdm/strategies/spectral.cpp
#         src/gdm/strategies/spectral_qawo.cpp
#         src/gdm/strategies/spectral_tanh_sinh.cpp
#         src/gdm/strategies/gorenflo_mainardi.cpp
#         src/gdm/strategies/huang_oberman.cpp
#         src/gdm/strategies/huang_oberman_linear.cpp
#         src/gdm/strategies/huang_oberman_quadratic.cpp
#         src/gdm/strategies/centered_3_point_periodized.cpp

# )


target_link_libraries(b118-frlap
    INTERFACE
        b118-core
        ${LIB_OPENBLAS}
        ${LIB_FFTW3}
)


# Debug Postffix
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(b118-frlap PROPERTIES DEBUG_POSTFIX "d")
endif()

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-frlap
    EXPORT "b118-frlap-targets"
    FILE_SET b118_frlap_public_headers
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"  # include
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" # include
)
