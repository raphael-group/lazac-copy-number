configure_file(
    "${lazac_SOURCE_DIR}/include/lazac.hpp.in"
    "${lazac_SOURCE_DIR}/include/lazac.hpp"
)

add_executable(lazac 
    lazac.cxx tree_io.cxx copy_number.cxx 
    clearcut.cxx cmdargs.cxx dist.cxx dmat.cxx fasta.cxx prng.cxx getopt_long.cxx
    )

# add libraries
target_link_libraries(lazac PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(lazac PRIVATE pprint)
target_link_libraries(lazac PRIVATE spdlog)
target_link_libraries(lazac PRIVATE csv)
target_link_libraries(lazac PRIVATE argparse)
