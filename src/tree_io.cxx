#include "tree_io.hpp"
#include <sstream>

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

    tree<float> read_newick_node(const std::string &newick, size_t &position) {
        tree<float> root(0);

        token t = read_token(newick, position);
        if (t != (token) separator::LEFT_PAREN){ // we are at a leaf
            t = read_token(newick, position);
            if (!std::holds_alternative<std::string>(t)) {
                root.name = "";
                position++;
            } else {
                root.name = std::get<std::string>(t);
                position += root.name.length();
            }

            return root;
        }

        position++;

        // otherwise, read until closing parentheses
        std::vector<tree<float>> children;
        while (position < newick.length()) {
            t = read_token(newick, position);
            if (t == (token) separator::LEFT_PAREN) {
                tree<float> child = read_newick_node(newick, position);
                children.push_back(child);
            } else if (t == (token) separator::COMMA) {
                position++;
                continue;
            } else if (t == (token) separator::RIGHT_PAREN){
                position++;
                break;
            } else {
                // issue is that stream is already partially consumed
                tree<float> child = read_newick_node(newick, position);
                children.push_back(child);
            }
        }

        for (const auto& child : children) {
            root.add_child(child);
        }

        t = read_token(newick, position);
        if (!std::holds_alternative<std::string>(t)) {
            position++;
            root.name = "";
        } else {
            root.name = std::get<std::string>(t);
            position += root.name.length();
        }

        return root;
    }

    tree<float> read_newick_node(const std::string &newick) {
        size_t position = 0;
        return read_newick_node(newick, position);
    }

    std::string print_newick_tree(const tree<float> &T) {
        if (!T.has_children()) return T.name;

        std::string newick("(");
        int counter = 0;
        for (const auto& child : T.get_children()) {
            if(counter != 0) newick += ",";
            newick += print_newick_tree(child);
            counter++;
        }

        newick += ")";
        newick += T.name;
        return newick;
    }

};
