#ifndef _TREE_IO_H
#define _TREE_IO_H

#include <optional>
#include <iostream>
#include <variant>
#include <string>
#include <stdexcept>
#include <optional>

#include "digraph.hpp"

namespace treeio {
    class malformed_parse_exception : public std::runtime_error {
    public:
        malformed_parse_exception(std::string reason) : std::runtime_error("Malformed input - failed to parse.\nReason: " + reason) { }
    };

    enum separator {
        LEFT_PAREN,
        RIGHT_PAREN,
        COMMA,
        COLON,
        SEMICOLON
    };

    struct newick_vertex_data {
        std::string name;
        std::optional<float> in_branch_length;
    };

    typedef std::string name;
    typedef std::variant<separator, name> token;

    /*
      Reads the next token, that is, either a name, a
      distance or a separation token. Throws a malformed
      exception if the string is improperly formed.
    */
    token read_token(const std::string &newick, size_t position);

    /*
      Recursive descent Newick parser based on 
      the following grammar:

         Tree -> Subtree ";"
         Subtree -> Leaf | Internal
         Leaf -> Name
         Internal -> "(" BranchSet ")" Name
         BranchSet -> Branch | Branch "," BranchSet
         Branch -> Subtree Length
         Name -> empty | string
         Length -> empty | ":" number

    */
    digraph<newick_vertex_data> read_newick_node(const std::string &newick);

    /*
      Requires: 
        - type T has a .name attribute for printing.
     */
    template<class T>
    std::string print_newick_tree(const digraph<T> &tree, int root) {
        if (tree.out_degree(root) == 0) {
            std::string newick = tree[root].data.name;

            if (tree[root].data.in_branch_length.has_value()) {
                newick += ":" + std::to_string(tree[root].data.in_branch_length.value());
            }

            return newick;
        }

        std::string newick("(");
        int counter = 0;
        for (const auto& child: tree.successors(root)) {
            if(counter != 0) newick += ",";
            newick += print_newick_tree(tree, child);
            counter++;
        }

        newick += ")";
        newick += tree[root].data.name;
        if (tree[root].data.in_branch_length.has_value()) {
            newick += ":" + std::to_string(tree[root].data.in_branch_length.value());
        }
        return newick;
    }

    template<class T>
    std::string print_newick_tree(const digraph<T> &tree) {
        return print_newick_tree(tree, 0);
    }
};

#endif
