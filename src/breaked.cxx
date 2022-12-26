#include <spdlog/spdlog.h>
#include <argparse/argparse.hpp>

#include <stdexcept>
#include <iostream>
#include <istream>
#include <stack>

#include "breaked.hpp"
#include "phylogeny.hpp"
#include "tree_io.hpp"

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program(
        "breaked",
        std::to_string(BREAKED_VERSION_MAJOR) + "." + std::to_string(BREAKED_VERSION_MINOR)
    );

    // program.add_argument("seed_tree")
    // .help("Seed tree in Newick format.");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    std::string newick_string("(((AAAAA,B,(C,D)E)F))");
    auto t = treeio::read_newick_node(newick_string);
    std::cout << treeio::print_newick_tree(t) << std::endl;

    return 0;
}
