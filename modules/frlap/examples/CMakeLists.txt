add_executable(b118-frlap-example-huaob2-canonic
    example_huaob2_canonic.cpp
)

target_link_libraries(b118-frlap-example-huaob2-canonic
    PRIVATE
        b118-frlap-miniapp-gdm-calculator-simple
)
