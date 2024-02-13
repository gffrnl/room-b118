# Fractional Laplacian module: b118::frlap

set(B118_MODULE_FRLAP_HEADERS
    frlap.hpp
    frlap/gdm.hpp
)

include(prepend_path_to_files.cmake)

prepend_path_to_files(B118_MODULE_FRLAP_HEADERS "${B188_INCLUDE_DIR}/b118" updated_headers)
set(B118_MODULE_FRLAP_HEADERS ${updated_headers} PARENT_SCOPE)