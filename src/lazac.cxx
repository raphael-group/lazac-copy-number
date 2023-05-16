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
#include "lazac.hpp"
#include "tree_io.hpp"

using namespace std;
using namespace copynumber;

using json = nlohmann::json;

/*
  Reads copy number profiles from a CSV
  file representation of the profiles.
*/
std::map<std::string, copynumber_profile> read_cn_profiles(std::string cn_profile_file) {
    std::map<std::string, copynumber_profile> cn_profiles;
    csv::CSVReader reader(cn_profile_file);
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

    return cn_profiles;
}

void do_distance(argparse::ArgumentParser distance) {
    std::map<std::string, copynumber_profile> cn_profiles = read_cn_profiles(distance.get<std::string>("cn_profile"));
    std::map<std::string, breakpoint_profile> bp_profiles;
    for (const auto &[name, cn_profile] : cn_profiles) {
        auto bp_profile = convert_to_breakpoint_profile(cn_profile, 2);
        bp_profiles[name] = bp_profile;
    }

    std::vector<std::string> names; 
    for (const auto& [name, _] : bp_profiles) {
        names.push_back(name);
    }

    std::sort(names.begin(), names.end());

    // we use a simple representation of a matrix
    // since it is only used for this purpose
    // and resize space for the elements
    std::vector<std::vector<int>> distance_matrix;
    distance_matrix.resize(names.size());
    for (std::vector<int>::size_type i = 0; i < names.size(); i++) {
        distance_matrix[i].resize(names.size());
    }

    // fill in the elements of the matrix
    for (std::vector<int>::size_type i = 0; i < names.size(); i++) {
        if (i % 100 == 0) {
            spdlog::info("Built {} out of {} rows of distance matrix.", i, names.size());
        }

        for (std::vector<int>::size_type j = i; j < names.size(); j++) {
            int dist = breakpoint_magnitude(bp_profiles[names[i]] - bp_profiles[names[j]]);
            distance_matrix[i][j] = dist;
            distance_matrix[j][i] = dist;
        }
    }

    std::ofstream matrix_output(distance.get<std::string>("-o") + "_dist_matrix.csv", std::ios::out);
    for (std::vector<int>::size_type i = 0; i < names.size(); i++) {
        if (i != 0) matrix_output << ",";
        matrix_output << names[i];
    }
    matrix_output << std::endl;

    for (std::vector<int>::size_type i = 0; i < names.size(); i++) {
        matrix_output << names[i];
        for (std::vector<int>::size_type j = 0; j < names.size(); j++) {
            matrix_output << "," << distance_matrix[i][j];
        }
        matrix_output << std::endl;
    }
}

void do_nni(argparse::ArgumentParser nni) {
    std::stringstream printer_stream;
    pprint::PrettyPrinter printer(printer_stream);
    printer.compact(true);

    std::ifstream in(nni.get<std::string>("seed_tree"));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string seed_tree_newick = buffer.str();

    digraph<treeio::newick_vertex_data> t = treeio::read_newick_node(seed_tree_newick);

    std::map<std::string, copynumber_profile> cn_profiles = read_cn_profiles(nni.get<std::string>("cn_profile"));
    std::map<std::string, breakpoint_profile> bp_profiles;
    std::vector<genomic_bin> sorted_bins;
    for (const auto &[name, cn_profile] : cn_profiles) {
        auto bp_profile = convert_to_breakpoint_profile(cn_profile, 2);
        bp_profiles[name] = bp_profile;
        sorted_bins = bp_profile.bins;
    }

    /* Creates rectilinear vertex data */
    digraph<breakpoint_vertex_data> breakpoint_tree;
    digraph<rectilinear_vertex_data> rectilinear_tree;
    for (auto u : t.nodes()) {
        breakpoint_vertex_data d;
        rectilinear_vertex_data r;

        const auto& vertex_data = t[u];
        d.name = vertex_data.data.name;
        r.name = vertex_data.data.name;

        if (t.out_degree(u) == 0) {
            d.breakpoint_profile = bp_profiles[d.name].profile;
            r.start = bp_profiles[d.name].profile;
            r.end = bp_profiles[d.name].profile;
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
    for (; counter < nni.get<int>("-i"); iteration++) {
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
        std::uniform_real_distribution<double> aggression_distrib(0, nni.get<double>("-a"));
        stochastic_nni(candidate_tree, gen, aggression_distrib(gen));

        digraph<rectilinear_vertex_data> updated_tree = hill_climb(candidate_tree, gen, nni.get<bool>("-g"));
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

    auto final_tree = ancestral_labeling(candidate_trees[candidate_trees.size() - 1], 0, sorted_bins);
    std::string newick_string = treeio::print_newick_tree(final_tree);
    newick_string += ";";

    std::ofstream newick_output(nni.get<std::string>("-o") + "_tree.newick", std::ios::out);
    newick_output << newick_string;

    std::ofstream info_output(nni.get<std::string>("-o") + "_info.json", std::ios::out);
    info_output << progress_information.dump();
}

int main(int argc, char *argv[])
{
    auto console_logger = spdlog::stdout_color_mt("lazac");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");


    argparse::ArgumentParser program(
        "lazac",
        std::to_string(LAZAC_VERSION_MAJOR) + "." + std::to_string(LAZAC_VERSION_MINOR)
    );

    argparse::ArgumentParser distance(
        "distance"
    );

    distance.add_description("Computes a distance matrix on copy number profiles");

    distance.add_argument("cn_profile")
        .help("copy number profile in CSV format");

    distance.add_argument("-o", "--output")
        .help("prefix of the output files")
        .required();

    argparse::ArgumentParser nni(
        "nni"
    );

    nni.add_description("Infers a copy number tree using NNI operations");

    nni.add_argument("cn_profile")
        .help("copy number profile in CSV format");

    nni.add_argument("seed_tree")
        .help("seed tree in Newick format.");

    nni.add_argument("-o", "--output")
        .help("prefix of the output files")
        .required();

    nni.add_argument("-a", "--aggression")
        .help("aggression of stochastic perturbation in (0, infinity)")
        .default_value(1.0)
        .scan<'g', double>();

    nni.add_argument("-i", "--iterations")
        .help("number of iterations to perform without improvement before stopping")
        .default_value(100)
        .scan<'d', int>();

    nni.add_argument("-g", "--greedy")
        .help("use greedy hill climbing strategy as opposed to full NNI neighborhood exploration")
        .default_value(false)
        .implicit_value(true);

    program.add_subparser(nni);
    program.add_subparser(distance);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(nni)) {
            std::cerr << nni;
        } else if (program.is_subcommand_used(distance)) {
            std::cerr << distance;
        } else {
            std::cerr << program;
        }

        std::exit(1);
    }

    if (program.is_subcommand_used(nni)) {
        do_nni(nni);
    } else if (program.is_subcommand_used(distance)) {
        do_distance(distance);
    } else {
        std::cerr << program;
    }

    return 0;
}
