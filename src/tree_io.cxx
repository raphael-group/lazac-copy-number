#include "tree_io.hpp"
#include <sstream>
#include <iostream>

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
            case ':':
                if (!empty) goto exit_loop;
                return separator::COLON;
            default:
                token_stream << char(c);
                empty = false;
            }
        }
    exit_loop: ;

        return token_stream.str();
    }

    void read_newick_node(const std::string &newick, size_t &position,
                          digraph<std::string>& tree, int root) {
        // Case 1: Check if leaf
        token t = read_token(newick, position);
        if (t != (token) separator::LEFT_PAREN){ 
            if (!std::holds_alternative<std::string>(t)) {
                tree[root].data = "";
                position++;
            } else {
                tree[root].data = std::get<std::string>(t);
                position += tree[root].data.length();
            }

            return;
        }

        position++;

        // Case 2: At internal node
        std::vector<digraph<std::string>> children;
        while (true) {
            t = read_token(newick, position);
            if (t == (token) separator::LEFT_PAREN) {
                int v = tree.add_vertex("internal");
                tree.add_edge(root, v);
                read_newick_node(newick, position, tree, v);
            } else if (t == (token) separator::COMMA) {
                position++;
                continue;
            } else if (t == (token) separator::RIGHT_PAREN){
                position++;
                break;
            } else {
                int v = tree.add_vertex("internal");
                tree.add_edge(root, v);
                read_newick_node(newick, position, tree, v);
            }
        }

        if (t != (token) separator::RIGHT_PAREN) {
            throw malformed_parse_exception("Expected right parentheses.");
        }

        t = read_token(newick, position);
        if (!std::holds_alternative<std::string>(t)) {
            tree[root].data = "";
        } else {
            tree[root].data = std::get<std::string>(t);
            position += tree[root].data.length();
        }
    }

    digraph<std::string> read_newick_node(const std::string &newick) {
        digraph<std::string> t;
        int root = t.add_vertex("root");
        size_t position = 0;
        read_newick_node(newick, position, t, root);
        return t;
    }
};
