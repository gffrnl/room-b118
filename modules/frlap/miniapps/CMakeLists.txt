add_library(b118-frlap-miniapp-gdm-calculator-simple)

target_sources(b118-frlap-miniapp-gdm-calculator-simple
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/b118/frlap/gdm/calculators/simple.cpp
)

target_include_directories(b118-frlap-miniapp-gdm-calculator-simple
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(b118-frlap-miniapp-gdm-calculator-simple
    PUBLIC
        b118-frlap
)