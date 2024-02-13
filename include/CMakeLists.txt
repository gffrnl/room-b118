function(prepend_path_to_files FILE_LIST PATH_PREFIX OUTPUT_LIST)
    set(updated_list "")
    foreach(file ${${FILE_LIST}})
        list(APPEND updated_list "${PATH_PREFIX}/${file}")
    endforeach()
    set(${OUTPUT_LIST} ${updated_list} PARENT_SCOPE)
endfunction()

set(B118_MODULE_BASIC_HEADERS
    type_traits.hpp
    numbers.hpp
    almost_equal.hpp
    equal_within_ulps.hpp
    linspace.hpp
)

set(B118_MODULE_LINALG_HEADERS
    linalg.hpp
)

set(B118_MODULE_FRLAP_HEADERS
    frlap.hpp
    frlap/gdm.hpp
)

# # Iterate over each filename in B118_MODULE_FRLAP_HEADERS
# set(updated_headers "")
# foreach(header ${B118_MODULE_FRLAP_HEADERS})
#     # Create the new path by prepending B188_INCLUDE_DIR/b118 to each filename
#     set(new_header "${B188_INCLUDE_DIR}/b118/${header}")
#     # Replace the original filename with the new path
#     list(APPEND updated_headers ${new_header})
# endforeach()

# # Replace B118_MODULE_FRLAP_HEADERS with updated list
# set(B118_MODULE_FRLAP_HEADERS ${updated_headers} PARENT_SCOPE)

# Basic module
prepend_path_to_files(B118_MODULE_BASIC_HEADERS "${B188_INCLUDE_DIR}/b118" updated_headers)
set(B118_MODULE_BASIC_HEADERS ${updated_headers} PARENT_SCOPE)

# # Linear algebra module: b118::linalg
prepend_path_to_files(B118_MODULE_LINALG_HEADERS "${B188_INCLUDE_DIR}/b118" updated_headers)
set(B118_MODULE_LINALG_HEADERS ${updated_headers} PARENT_SCOPE)

# Fractional Laplacian module: b118::frlap
prepend_path_to_files(B118_MODULE_FRLAP_HEADERS "${B188_INCLUDE_DIR}/b118" updated_headers)
set(B118_MODULE_FRLAP_HEADERS ${updated_headers} PARENT_SCOPE)
