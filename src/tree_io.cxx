#include "tree_io.hpp"
#include <sstream>
#include <iostream>
#include <regex>

namespace treeio {
    token read_token(const std::string &newick, size_t position) {
        bool empty = true;
        std::stringstream token_stream;

        for(; position < newick.length(); position++) {
            char c = newick[position];
            
            switch(c) {
            case '(':
                if (!empty) goto exit_loop;
                return separator::LEFT_PAREN;
            case ')':
                if (!empty) goto exit_loop;
                return separator::RIGHT_PAREN;
            case ',':
                if (!empty) goto exit_loop;
                return separator::COMMA;
            case ';':
                if (!empty) goto exit_loop;
                return separator::SEMICOLON;
            default:
                token_stream << char(c);
                empty = false;
            }
        }
    exit_loop: ;

        return token_stream.str();
    }

    /*
       Parses a string of the form [^:]+:<float> into a tuple (name, length).
     */
    newick_vertex_data parse_newick_vertex(std::string s) {
        std::regex rgx("([^:]*):?(.*)");
        std::smatch match;

        if (!std::regex_match(s, match, rgx)) {
            throw malformed_parse_exception("Malformed vertex name.");
        }

        newick_vertex_data d;
        d.name = match[1];
        if (match[2] != "") {
            d.in_branch_length = std::stof(match[2]);
        }

        return d;
    }

    void read_newick_node(const std::string &newick, size_t &position,
                          digraph<newick_vertex_data>& tree, int root, int &internal_counter) {
        // Case 1: Check if leaf
        token t = read_token(newick, position);
        if (t != (token) separator::LEFT_PAREN){ 
            if (!std::holds_alternative<std::string>(t)) {
                position++;
            } else {
                std::string vertex_string = std::get<std::string>(t);
                tree[root].data = parse_newick_vertex(vertex_string);
                position += vertex_string.length();
            }

            return;
        }

        position++;

        // Case 2: At internal node
        std::vector<digraph<std::string>> children;
        while (true) {
            t = read_token(newick, position);
            if (t == (token) separator::LEFT_PAREN) {
                newick_vertex_data d;
                d.name = "internal_" + std::to_string(internal_counter++);
                int v = tree.add_vertex(d);
                tree.add_edge(root, v);
                read_newick_node(newick, position, tree, v, internal_counter);
            } else if (t == (token) separator::COMMA) {
                position++;
                continue;
            } else if (t == (token) separator::RIGHT_PAREN){
                position++;
                break;
            } else {
                newick_vertex_data d;
                d.name = "";
                int v = tree.add_vertex(d);
                tree.add_edge(root, v);
                read_newick_node(newick, position, tree, v, internal_counter);
            }
        }

        if (t != (token) separator::RIGHT_PAREN) {
            throw malformed_parse_exception("Expected right parentheses.");
        }

        t = read_token(newick, position);
        if (std::holds_alternative<std::string>(t)) {
            std::string vertex_string = std::get<std::string>(t);
            tree[root].data = parse_newick_vertex(vertex_string);
            position += vertex_string.length();
        }
    }

    digraph<newick_vertex_data> read_newick_node(const std::string &newick) {
        digraph<newick_vertex_data> t;
        newick_vertex_data d;
        d.name = "root";
        int root = t.add_vertex(d);
        size_t position = 0;
        int counter = 0;
        read_newick_node(newick, position, t, root, counter);
        return t;
    }

};
