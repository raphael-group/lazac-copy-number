configure_file(
    "${breaked_SOURCE_DIR}/include/breaked.h.in" 
    "${breaked_SOURCE_DIR}/include/breaked.h"
)

add_executable(breaked breaked.cxx)

# add libraries
target_link_libraries(breaked PRIVATE spdlog)
target_link_libraries(breaked PRIVATE csv)
target_link_libraries(breaked PRIVATE argparse)
