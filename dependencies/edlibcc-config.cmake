add_library(EdlibCC STATIC
        ${CMAKE_CURRENT_LIST_DIR}/edlib/edlib.cc)

set_target_properties(EdlibCC PROPERTIES
    INTERFACE_COMPILE_OPTIONS "-std=c++11"
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/edlib/"
)