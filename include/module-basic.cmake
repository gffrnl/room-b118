# Basic module

add_library(b118-basic)

target_sources(b118-basic
    PUBLIC
        FILE_SET b118_basic_public_headers
            TYPE HEADERS
            FILES
                b118/type_traits.hpp
                b118/numbers.hpp
                b118/almost_equal.hpp
                b118/equal_within_ulps.hpp
                b118/linspace.hpp
                # vendors
                b118/vendor/cppreference/equal_within_ulps.hpp
)

target_sources(b118-basic PRIVATE src/file.cpp)

target_include_directories(b118-basic PUBLIC "${B188_INCLUDE_DIR}")

# Debug Suffix
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(b118-basic PROPERTIES DEBUG_POSTFIX "d")
endif()

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-basic
  EXPORT "b118-basic-targets"
  FILE_SET b118_basic_public_headers
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)
