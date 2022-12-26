#ifndef _TREE_IO_H
#define _TREE_IO_H

#include <variant>
#include <string>
#include <stdexcept>

#include "phylogeny.hpp"

namespace treeio {
    class malformed_parse_exception : public std::runtime_error {
    public:
        malformed_parse_exception() : std::runtime_error("Malformed input - failed to parse.") { }
    };

    enum separator {
        LEFT_PAREN,
        RIGHT_PAREN,
        COMMA,
        COLON,
        SEMICOLON
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
    tree<float> read_newick_node(const std::string &newick, size_t &position);
    tree<float> read_newick_node(const std::string &newick);

    std::string print_newick_tree(const tree<float> &T);
};

#endif
