
add_library(b118-linalg)

target_sources(b118-linalg
    PUBLIC
        FILE_SET b118_linalg_public_headers
            TYPE HEADERS
            FILES
                b118/linalg.hpp
)

target_sources(b118-linalg
    PRIVATE
        src/b118-linalg.cpp
)

target_include_directories(b118-linalg
    PUBLIC
        ${CMAKE_SOURCE_DIR}/modules/base
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(b118-linalg
    PUBLIC
        b118-base
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
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)
