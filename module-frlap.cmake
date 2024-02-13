# Fractional Laplacian module: b118::frlap

add_library(b118-frlap)

target_sources(b118-frlap
    PUBLIC
        FILE_SET b118_frlap_public_headers
            TYPE HEADERS
            FILES
                ${B188_INCLUDE_DIR}/b118/linspace.hpp
)

target_sources(b118-frlap PRIVATE src/file.cpp)

target_include_directories(b118-frlap PUBLIC "${B188_INCLUDE_DIR}")

# Debug Suffix
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
