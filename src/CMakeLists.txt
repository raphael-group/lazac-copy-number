configure_file(
    "${breaked_SOURCE_DIR}/include/breaked.hpp.in" 
    "${breaked_SOURCE_DIR}/include/breaked.hpp"
)

add_executable(breaked breaked.cxx tree_io.cxx copy_number.cxx)

# add libraries
target_link_libraries(breaked PRIVATE pprint)
target_link_libraries(breaked PRIVATE spdlog)
target_link_libraries(breaked PRIVATE csv)
target_link_libraries(breaked PRIVATE argparse)
