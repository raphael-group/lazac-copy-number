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
        succ[u].insert(v);
        pred[v].insert(u);
    }

    void remove_edge(int u, int v) {
        succ[u].erase(v);
        pred[v].erase(u);
    }

    std::vector<int> nodes() const {
        std::vector<int> vertices;
        for (int i = 0; i < id_counter; i++) {
            vertices.push_back(i);
        }
        return vertices;
    }

    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> edges;
        for (const auto &[u, vs] : succ) {
            for (const auto &v : vs) {
                edges.push_back(std::make_pair(u, v));
            }
        }
        return edges;
    }

    /* TODO
       WARNING: does not maintain invariant
       that all edges are between 0 and N. We 
       will need to update the code to fix this.
       Probably should return a map from old to new
       vertex labels.
     */
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

    const std::set<int>& predecessors(int u) const {
        return pred.at(u);
    }

    const std::set<int>& successors(int u) const {
        return succ.at(u);
    }

    size_t out_degree(int u) const {
        return succ.at(u).size();
    }
};

#endif
