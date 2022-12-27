#ifndef _DIGRAPH_H
#define _DIGRAPH_H

#include <vector>
#include <set>
#include <map>
#include <iostream>

/*
  Simple adjacency list representation of
  a DiGraph.
*/

template <class T>
class vertex {
public:
    int id;
    T data;
    vertex() {};
    vertex(int id, T data) : id(id), data(data) {};
};

template <class T>
class digraph {
private:
    int id_counter = 0;

    std::map<int, std::set<int>> succ;
    std::map<int, std::set<int>> pred;

    std::map<int, vertex<T>> vertices;
public:
    // returns id of created vertex
    int add_vertex(T data) {
        vertex<T> v(id_counter, data);
        vertices[v.id] = v;
        succ[v.id] = std::set<int>();
        pred[v.id] = std::set<int>();
        id_counter++;
        return v.id;
    }

    void add_edge(int u, int v) {
        // assert u and v in graph
        succ[u].insert(v);
        pred[v].insert(u);
    }

    void delete_vertex(int u) {
        // removes u and all (v, u) edges
        vertices.erase(u);
        for (const vertex<T>& v : pred[u]) {
            succ[v].erase(u);
        }

        succ.erase(u);
        pred.erase(u);
    }

    vertex<T>& operator[](int u) {
        return vertices.at(u);
    }

    const vertex<T>& operator[](int u) const {
        return vertices.at(u);
    }

    const std::set<int>& neighbors(int u) const {
        return succ.at(u);
    }

    size_t out_degree(int u) const {
        return succ.at(u).size();
    }

    //std::string to_string() const {
    //std::string out;
    //for (const auto& kv : succ) {
    //out += kv.
        //}
    //
    //}
};

#endif

