#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include <random>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <pprint.hpp>
#include <optional>
#include <stack>
#include <tuple>

#include "copy_number.hpp"
#include "digraph.hpp"
#include "breaked.hpp"
#include "tree_io.hpp"

using namespace std;
using namespace copynumber;

int main(int argc, char *argv[])
{
    auto console_logger = spdlog::stdout_color_mt("breaked");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");

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

    pprint::PrettyPrinter printer;
    std::map<std::string, breakpoint_profile> bp_profiles;
    for (const auto &[name, cn_profile] : cn_profiles) {
        auto bp_profile = convert_to_breakpoint_profile(cn_profile, 2);
        bp_profiles[name] = bp_profile;
    }

    /* Creates rectilinear vertex data */
    digraph<breakpoint_vertex_data> breakpoint_tree;
    digraph<rectilinear_vertex_data> rectilinear_tree;
    for (auto u : t.nodes()) {
        breakpoint_vertex_data d;
        rectilinear_vertex_data r;

        const auto &name = t[u];
        d.name = name.data;
        r.name = name.data;

        if (name.data != "") {
            d.breakpoint_profile = bp_profiles[name.data].profile;
            r.start = bp_profiles[name.data].profile;
            r.end = bp_profiles[name.data].profile;
        } 

        breakpoint_tree.add_vertex(d);
        rectilinear_tree.add_vertex(r);
    }

    for (auto [u, v] : t.edges()) {
        breakpoint_tree.add_edge(u, v);
        rectilinear_tree.add_edge(u, v);
    }

    std::random_device rd;
    std::ranlux48_base gen(rd());
    for (int i = 0; i < 2; i++) {
        small_rectilinear(rectilinear_tree, 0);
        std::cout << rectilinear_tree[0].data.score << std::endl;
        unvisit(rectilinear_tree, 0);

        digraph<rectilinear_vertex_data> t2 = stochastic_nni(rectilinear_tree, gen, 0.25);
        std::cout << "here" << std::endl;
        small_rectilinear(t2, 0);
        std::cout << t2[0].data.score << std::endl;
    }

    return 0;
}
