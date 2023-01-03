#include <pprint.hpp>
#include <nlohmann/json.hpp>
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
#include <optional>
#include <stack>
#include <tuple>

#include "copy_number.hpp"
#include "digraph.hpp"
#include "breaked.hpp"
#include "tree_io.hpp"

using namespace std;
using namespace copynumber;

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    auto console_logger = spdlog::stdout_color_mt("breaked");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");

    std::stringstream printer_stream;
    pprint::PrettyPrinter printer(printer_stream);
    printer.compact(true);

    argparse::ArgumentParser program(
        "breaked",
        std::to_string(BREAKED_VERSION_MAJOR) + "." + std::to_string(BREAKED_VERSION_MINOR)
    );

    program.add_argument("cn_profile")
        .help("copy number profile in CSV format");

    program.add_argument("seed_tree")
        .help("seed tree in Newick format.");

    program.add_argument("-o", "--output")
        .help("prefix of the output files")
        .required();

    program.add_argument("-a", "--aggression")
        .help("aggression of stochastic perturbation in (0, infinity)")
        .default_value(1.0)
        .scan<'g', double>();

    program.add_argument("-i", "--iterations")
        .help("number of iterations to perform without improvement before stopping")
        .default_value(100)
        .scan<'d', int>();

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

    /*
      Candidate tree set is obtained by randomly
      perturbing candidate trees.
     */
    std::vector<digraph<rectilinear_vertex_data>> candidate_trees;
    float aggressions[] = {0, 0.25, 0.50, 0.75, 1, 1.25, 1.5, 1.75};
    for (float aggression : aggressions) {
        digraph<rectilinear_vertex_data> t = stochastic_nni(rectilinear_tree, gen, aggression);
        candidate_trees.push_back(t);
    }

    json progress_information;
    int counter = 0, iteration = 0;
    for (; counter < program.get<int>("-i"); iteration++) {
        for (auto& candidate_tree : candidate_trees) {
            small_rectilinear(candidate_tree, 0);
        }

        std::sort(candidate_trees.begin(), candidate_trees.end(),
                  [](const digraph<rectilinear_vertex_data> &a, const digraph<rectilinear_vertex_data> &b) {
                      return a[0].data.score > b[0].data.score;
        });

        std::vector<int> scores;
        for (auto& candidate_tree : candidate_trees) {
            scores.push_back(candidate_tree[0].data.score);
        }

        json progress_information_i;
        progress_information_i["scores"] = scores;
        progress_information_i["iteration"] = iteration;
        progress_information.push_back(progress_information_i);

        printer.print(scores);
        std::string candidate_scores_string = printer_stream.str();
        candidate_scores_string = candidate_scores_string.substr(0, candidate_scores_string.length() - 1);
        spdlog::info("Candidate tree scores @ iteration {}: {}", iteration, candidate_scores_string);
        printer_stream.str("");

        // Select and perturb candidate tree.
        std::uniform_int_distribution<int> distrib(0, candidate_trees.size() - 1);
        int candidate_tree_idx = distrib(gen);

        digraph<rectilinear_vertex_data> candidate_tree = candidate_trees[candidate_tree_idx];
        stochastic_nni(candidate_tree, gen, program.get<double>("-a"));

        digraph<rectilinear_vertex_data> updated_tree = hill_climb(candidate_tree);
        if (updated_tree[0].data.score < candidate_trees[0][0].data.score) {
            candidate_trees[0] = updated_tree;
            spdlog::info("Updated candidate tree set.");
            counter = 0;
            continue;
        } 

        counter++;
    }

    for (auto& candidate_tree : candidate_trees) {
        small_rectilinear(candidate_tree, 0);
    }

    std::sort(candidate_trees.begin(), candidate_trees.end(),
              [](const digraph<rectilinear_vertex_data> &a, const digraph<rectilinear_vertex_data> &b) {
                  return a[0].data.score > b[0].data.score;
              });

    std::string newick_string = treeio::print_newick_tree(candidate_trees[candidate_trees.size() - 1]);
    newick_string += ";";

    std::ofstream newick_output(program.get<std::string>("-o") + "_tree.newick", std::ios::out);
    newick_output << newick_string;

    std::ofstream info_output(program.get<std::string>("-o") + "_info.json", std::ios::out);
    info_output << progress_information.dump();

    return 0;
}
