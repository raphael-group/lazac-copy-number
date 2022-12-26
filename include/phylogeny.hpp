#ifndef _PHYLOGENY_H
#define _PHYLOGENY_H

#include <string>
#include <vector>

template <class T>
class tree {
private:
    std::vector<tree<T>> children;
    tree<T> *parent = nullptr;

public:
    std::string name;
    T data;

    tree(T data, std::string name) : data(data), name(name) {};
    tree(T data) : data(data), name("") {};

    void add_child(tree<T> child) {
        child.set_parent(this);
        this->children.push_back(child);
    }

    void set_parent(tree<T> *parent) {
        this->parent = parent;
    }

    tree<T>* get_parent() const {
        return this->parent;
    }

    const std::vector<tree<T>>& get_children() const {
        return this->children;
    }

    bool has_children() const {
        return !this->children.empty();
    }
};

#endif
