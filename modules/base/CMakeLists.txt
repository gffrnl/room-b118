
add_library(b118-base)

target_sources(b118-base
    PUBLIC
        FILE_SET b118_base_public_headers
            TYPE HEADERS
            FILES
                b118/utility.hpp
                b118/type_traits.hpp
                b118/numbers.hpp
                b118/almost_equal.hpp
                b118/equal_within_ulps.hpp
                b118/closest_sorted.hpp
                b118/linspace.hpp
                # vendors
                b118/vendor/cppreference/equal_within_ulps.hpp
)

target_sources(b118-base
    PRIVATE
        src/b118-base.cpp
)

target_include_directories(b118-base
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_definitions(b118-base
    PUBLIC B118_VERSION="${CMAKE_PROJECT_VERSION}"
    PUBLIC B118_VERSION_MAJOR="${CMAKE_PROJECT_VERSION_MAJOR}"
    PUBLIC B118_VERSION_MINOR="${CMAKE_PROJECT_VERSION_MINOR}"
    PUBLIC B118_VERSION_PATCH="${CMAKE_PROJECT_VERSION_PATCH}"
    PUBLIC B118_VERSION_TWEAK="${CMAKE_PROJECT_VERSION_TWEAK}"
)

# Debug Postffix
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(b118-base PROPERTIES DEBUG_POSTFIX "d")
endif()

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-base
  EXPORT "b118-base-targets"
  FILE_SET b118_base_public_headers
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)

#add_subdirectory(examples)
