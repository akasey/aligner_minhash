add_library(SswCC::Static INTERFACE IMPORTED)

set_target_properties(SswCC::Static PROPERTIES
        INTERFACE_COMPILE_OPTIONS "-std=c11"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_LIST_DIR}/ssw/src/"
        )