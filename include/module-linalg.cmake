# Linear algebra module: b118::linalg

set(B118_MODULE_LINALG_HEADERS
    linalg.hpp
)

include(prepend_path_to_files.cmake)

prepend_path_to_files(B118_MODULE_LINALG_HEADERS "${B188_INCLUDE_DIR}/b118" updated_headers)
set(B118_MODULE_LINALG_HEADERS ${updated_headers} PARENT_SCOPE)