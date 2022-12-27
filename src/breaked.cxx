#include <spdlog/spdlog.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include <fstream>
#include <stdexcept>
#include <iostream>

#include "digraph.hpp"
#include "breaked.hpp"
#include "phylogeny.hpp"
#include "tree_io.hpp"

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program(
        "breaked",
        std::to_string(BREAKED_VERSION_MAJOR) + "." + std::to_string(BREAKED_VERSION_MINOR)
    );

    program.add_argument("seed_tree")
        .help("Seed tree in Newick format.");

    program.add_argument("cn_profile")
        .help("Copy number profile in CSV format.");

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    std::ifstream in(program.get<std::string>("seed_tree"));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string seed_tree_newick = buffer.str();

    digraph<std::string> t = treeio::read_newick_node(seed_tree_newick);
    std::cout << treeio::print_newick_tree(t) << std::endl;

    csv::CSVReader reader(program.get<std::string>("cn_profile"));
    for (const csv::CSVRow &row: reader) {
        // std::cout << row["node"].get<std::string>() << std::endl;
        // std::cout << row["chrom"].get<std::string>() << std::endl;
    }

    return 0;
}
