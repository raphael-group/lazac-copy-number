#include <spdlog/spdlog.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <pprint.hpp>

#include "copy_number.hpp"
#include "digraph.hpp"
#include "breaked.hpp"
#include "tree_io.hpp"

using namespace std;
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

    std::map<std::string, copynumber_profile> cn_profiles;
    csv::CSVReader reader(program.get<std::string>("cn_profile"));
    for (const csv::CSVRow &row: reader) {
        std::string node = row["node"].get<std::string>();
        std::string chrom = row["chrom"].get<std::string>();
        int start = row["start"].get<int>();
        int end = row["end"].get<int>();
        int cn_a = row["cn_a"].get<int>();

        genomic_bin bin(chrom, "cn_a", start, end);
        if (!cn_profiles.count(node)) {
            cn_profiles[node] = copynumber_profile();
        }

        cn_profiles[node].profile.push_back(cn_a);
        cn_profiles[node].bins.push_back(bin);
    }

    //pprint::PrettyPrinter printer;
    std::map<std::string, breakpoint_profile> bp_profiles;
    for (auto &[name, cn_profile] : cn_profiles) {
        auto bp_profile = convert_to_breakpoint_profile(cn_profile, 2);
        bp_profiles[name] = bp_profile;
    }

    return 0;
}
