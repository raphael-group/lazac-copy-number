#ifndef _COPY_NUMBER_H
#define _COPY_NUMBER_H

#include "digraph.hpp"

#include <random>
#include <vector>
#include <string>
#include <optional>
#include <stack>
#include <tuple>

namespace copynumber {
    struct genomic_bin {
        std::string chromosome;
        std::string allele;
        int start;
        int end;

        genomic_bin() {};
        genomic_bin(std::string chromosome, std::string allele, int start, int end) :
            chromosome(chromosome), allele(allele), start(start), end(end) {};

        bool operator==(const genomic_bin& other) const {
            return (chromosome == other.chromosome &&
                    allele == other.allele &&
                    start == other.start &&
                    end == other.end);
        }

        /* Lexicographic ordering */
        bool operator<(const genomic_bin& other) const {
            if (chromosome < other.chromosome) return true;
            if (chromosome > other.chromosome) return false;
            if (allele < other.allele) return true;
            if (allele > other.allele) return false;
            if (start < other.start) return true;
            if (start > other.start) return false;
            if (end < other.end) return true;
            if (end > other.end) return false;
            return false;
        }

        friend std::ostream& operator<<(std::ostream& os, const genomic_bin& bin);
    };

    struct copynumber_profile {
        std::vector<int> profile;
        std::vector<genomic_bin> bins;

        copynumber_profile() {};
        copynumber_profile(std::vector<int> profile, std::vector<genomic_bin> bins) :
            profile(profile), bins(bins) {};
    };

    struct breakpoint_profile {
        std::vector<int> profile;
        std::vector<genomic_bin> bins;

        breakpoint_profile() {};
        breakpoint_profile(std::vector<int> profile, std::vector<genomic_bin> bins) :
            profile(profile), bins(bins) {};

        breakpoint_profile operator+(const breakpoint_profile& other) const {
            breakpoint_profile result;
            result.profile = std::vector<int>(this->profile.size());
            for (std::vector<int>::size_type i = 0; i < this->profile.size(); i++) {
                result.profile[i] = this->profile[i] + other.profile[i];
            }
            result.bins = this->bins;
            return result;
        }

        breakpoint_profile operator-(const breakpoint_profile& other) const {
            breakpoint_profile result;
            result.profile = std::vector<int>(this->profile.size());
            for (std::vector<int>::size_type i = 0; i < this->profile.size(); i++) {
                result.profile[i] = this->profile[i] - other.profile[i];
            }
 
            result.bins = this->bins;
            return result;
        }
    };

    struct breakpoint_vertex_data {
        std::string name;
        std::optional<std::vector<int>> breakpoint_profile;
    };

    /*
      Rectilinear invariant: If visited == true for some vertex u,
      then for all i, the interval [start_i, end_i] is the set of
      values minimizing the rectilinear score of the sub-tree rooted
      at u and score is the minimizing rectlinear score for that
      sub-tree.
     */
    struct rectilinear_vertex_data {
        std::string name;
        std::optional<std::vector<int>> start;
        std::optional<std::vector<int>> end;

        int score = 0;
        bool visited = false;
    };

    /*
      Creates a breakpoint profile from a copy number profile
      where each chromosome and allele pair is considered seperately
      and genomic bins are sorted by their starting position.
    */
    breakpoint_profile convert_to_breakpoint_profile(const copynumber_profile &p, int diploid_cn);

    /*
      Overlaps two intervals [s1, e1], [s2, e2], returning the empty set 
      if they do not overlap.
     */
    std::optional<std::pair<int, int>> overlap(int s1, int e1, int s2, int e2);
    std::tuple<std::vector<int>, std::vector<int>, int> sankoff(const rectilinear_vertex_data& u, const rectilinear_vertex_data& v);

    /*
      Solves the small rectilinear problem for the sub-trees
      rooted at every vertex. Avoids excess recomputation.
      
      Requires:
        - t satisfies the *rectilinear invariant*.
         
      Output guarantees:
        - sets visited == true for all vertices in t.
        - t satisfies the *rectilinear invariant*.
    */
    void small_rectilinear(digraph<rectilinear_vertex_data>& t, int root);

    /*
      Performs (or undos) a NNI operation on edges (u, w) and (v, z) by
      swapping the edges.
      
      Does not necessarily maintain rectlinear invariant, but it can be
      re-maintained by calling unvisit(t, u) and unvisit(t, w).
    */
    void nni(digraph<rectilinear_vertex_data>& t, int u, int w, int v, int z);
    void undo_nni(digraph<rectilinear_vertex_data>& t, int u, int w, int v, int z);

    /*
      Unvisits all vertices on path from root to u. Trivially
      guarantees the *rectilinear invariant*.
     */
    void unvisit(digraph<rectilinear_vertex_data> &t, int root, int u);

    /*
      Unvisits all vertices in sub-tree rooted at root. Trivially
      guarantees the *rectilinear invariant*.
     */
    void unvisit(digraph<rectilinear_vertex_data> &t, int root);

    /*
      Assumes the root is the 0 vertex.
    */
    digraph<rectilinear_vertex_data> stochastic_nni(const digraph<rectilinear_vertex_data>& t, std::ranlux48_base& gen, float aggression);


    /*
    * Performs hill climbing on the rectilinear score of the input tree until no
    * more improvement is found.
    *
    * Parameters
    *    - t: input tree does not need to satisfy rectilinear invariant.
    *    - greedy: if true, selects the first improvement at every iteration. otherwise, 
    *      explores entire NNI neighborhood for improvement at every iteration.
    *    - gen: random generator to shuffle edges for random exploration.
    */
    digraph<rectilinear_vertex_data> hill_climb(digraph<rectilinear_vertex_data> t, std::ranlux48_base& gen, bool greedy);

    /*
      Computes the breakpoint magnitude of a *chromosome and allele sorted*
      chromosome breakpoint profile.
     */
    int breakpoint_magnitude(const breakpoint_profile& p);
};

#endif
