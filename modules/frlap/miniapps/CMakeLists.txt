add_library(b118-frlap-miniapp-gdm-calculator-simple)

target_sources(b118-frlap-miniapp-gdm-calculator-simple
    PUBLIC
        FILE_SET b118_frlap_miniapp_gdm_calculator_simple_public_headers
            TYPE HEADERS
            FILES
                b118/frlap/gdm/calculators/simple.hpp
)


target_sources(b118-frlap-miniapp-gdm-calculator-simple
    PRIVATE
        b118/frlap/gdm/calculators/simple.cpp
)

target_include_directories(b118-frlap-miniapp-gdm-calculator-simple
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(b118-frlap-miniapp-gdm-calculator-simple
    PUBLIC
        b118-frlap
)

# Debug Postffix
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_target_properties(b118-frlap PROPERTIES DEBUG_POSTFIX "d")
endif()

# for CMAKE_INSTALL_INCLUDEDIR definition
include(GNUInstallDirs)

# install the target and create export-set
install(TARGETS b118-frlap-miniapp-gdm-calculator-simple
  EXPORT "b118-frlap-miniapp-gdm-calculator-simple-targets"
  FILE_SET b118_frlap_miniapp_gdm_calculator_simple_public_headers
  # these get default values from GNUInstallDirs, no need to set them
  #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
  #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
  # except for public headers, as we want them to be inside a library folder
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)
