find_library(LIB_OPENBLAS NAMES openblas REQUIRED)
find_library(LIB_FFTW3 NAMES fftw3 REQUIRED)

message(STATUS "LIB_OPENBLAS = ${LIB_OPENBLAS}")
message(STATUS "LIB_FFTW3 = ${LIB_FFTW3}")

add_library(b118-frlap)

target_sources(b118-frlap
    PUBLIC
        FILE_SET b118_frlap_public_headers
            TYPE HEADERS
            FILES
                b118/frlap.hpp
                b118/frlap/gdm.hpp
                b118/frlap/gdm/strategy.hpp
                b118/frlap/gdm/trunc_uniform.hpp
                # strategies
                b118/frlap/gdm/strategies/spectral_qawo.hpp
                b118/frlap/gdm/strategies/spectral_tanh_sinh.hpp
                b118/frlap/gdm/strategies/gorenflo_mainardi.hpp
                b118/frlap/gdm/strategies/huang_oberman_linear.hpp
                b118/frlap/gdm/strategies/huang_oberman_quadratic.hpp
                b118/frlap/gdm/strategies/centered_3_point_periodized.hpp
                # toep.h (TO REMOVE)
                b118/frlap/detail/toep.h # TODO(gffrnl) to remove
                b118/frlap/detail/fstp.hpp
                b118/frlap/gdm/convolucao.hpp
)

target_sources(b118-frlap
    PRIVATE
        src/b118-frlap.cpp
        # strategies
        src/gdm/strategies/spectral_qawo.cpp
        src/gdm/strategies/spectral_tanh_sinh.cpp
        src/gdm/strategies/gorenflo_mainardi.cpp
        src/gdm/strategies/huang_oberman.cpp
        src/gdm/strategies/huang_oberman_linear.cpp
        src/gdm/strategies/huang_oberman_quadratic.cpp
        src/gdm/strategies/centered_3_point_periodized.cpp
    PRIVATE
        #src/toep.c
        src/fstp.cpp
)

target_include_directories(b118-frlap
    PUBLIC
        ${CMAKE_SOURCE_DIR}/modules/base
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(b118-frlap
    PUBLIC
        b118-base
    PUBLIC
        ${LIB_OPENBLAS}
        ${LIB_FFTW3}
)

add_dependencies(b118-frlap b118-base)

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
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)

add_subdirectory(miniapps)
#add_subdirectory(examples)
