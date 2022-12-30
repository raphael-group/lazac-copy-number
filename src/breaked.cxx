#include <spdlog/spdlog.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include <tuple>
#include <stack>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <pprint.hpp>
#include <optional>

#include "copy_number.hpp"
#include "digraph.hpp"
#include "breaked.hpp"
#include "tree_io.hpp"

struct breakpoint_vertex_data {
    std::string name;
    std::optional<std::vector<int>> breakpoint_profile;
};

struct rectilinear_vertex_data {
    std::string name;
    std::optional<std::vector<int>> start;
    std::optional<std::vector<int>> end;

    int score = 0;
    bool visited = false;
};

std::optional<std::pair<int, int>> overlap(int s1, int e1, int s2, int e2) {
    std::optional<std::pair<int, int>> out_interval;
    // 2 == self && 1 == other
    if (s1 <= e2 && s1 >= s2) {
        out_interval = std::make_pair(s1, std::min(e1, e2));
    } else if (s2 <= e1 and s2 >= s1) {
        out_interval = std::make_pair(s2, std::min(e1, e2));
    }

    return out_interval;
}

std::tuple<std::vector<int>, std::vector<int>, int> sankoff(const rectilinear_vertex_data& u, const rectilinear_vertex_data& v) {
    const std::vector<int> &u_start = u.start.value();
    const std::vector<int> &u_end = u.end.value();
    const std::vector<int> &v_start = v.start.value();
    const std::vector<int> &v_end = v.end.value();

    std::vector<int> start(u_start.size());
    std::vector<int> end(u_end.size());
    int distance = 0;
    for (int i = 0; i < u_start.size(); i++) {
        std::optional<std::pair<int, int>> interval = overlap(u_start[i], u_end[i], v_start[i], v_end[i]);

        if (interval) {
            auto [s, e] = interval.value();
            start[i] = s;
            end[i] = e;
        } else {
            if (u_start[i] < v_start[i]) {
                start[i] = u_end[i];
                end[i] = v_start[i];
            } else {
                start[i] = v_end[i];
                end[i] = u_start[i];
            }

            distance += end[i] - start[i];
        }
    }

    return std::make_tuple(start, end, distance);
}

/*
  Modifies the internal nodes of t but does not modify the
  leaves.
 */
void small_rectilinear(digraph<rectilinear_vertex_data>& t, int root) {
    std::stack<int> callstack;

    callstack.push(root);
    while (!callstack.empty()) {
        int node = callstack.top();
        callstack.pop();

        if (t.out_degree(node) == 0) {
            t[node].data.visited = true;
            continue;
        }

        // check condition that every node has two children
        if (t.out_degree(node) != 2)
            throw std::logic_error("every child must have exactly two children");

        // check to see if all children are visited
        bool children_visited = true; 
        for (const auto& child : t.neighbors(node)) {
            if (!t[child].data.visited) children_visited = false;
        }

        if (children_visited) {
            // grab two children
            std::set<int> children_s = t.neighbors(node);
            std::vector<int> children_a(children_s.begin(), children_s.end());

            /* could wrap this into "Sankoff" sub-routine */
            int u = children_a[0];
            int v = children_a[1];

            const rectilinear_vertex_data& u_data = t[u].data;
            const rectilinear_vertex_data& v_data = t[v].data;

            const auto& [start, end, cost] = sankoff(u_data, v_data);

            t[node].data.score = cost + u_data.score + v_data.score;
            t[node].data.start = start;
            t[node].data.end = end;
            t[node].data.visited = true;

            continue;
        }

        callstack.push(node);
        for (const auto& child : t.neighbors(node)) {
            callstack.push(child);
        }
    }
}

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

    small_rectilinear(rectilinear_tree, 0);
    for (const auto& child : rectilinear_tree.neighbors(0)) {
    }

    return 0;
}
