#include <spdlog/spdlog.h>
#include <argparse/argparse.hpp>

#include <iostream>

#include "breaked.h"

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("breaked");

    program.add_argument("seed_tree")
        .help("Seed tree")
        .scan<'i', int>();

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    return 0;
}
