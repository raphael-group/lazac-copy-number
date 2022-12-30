#ifndef _VEC_UTILS_H
#define _VEC_UTILS_H

#include <numeric>
#include <vector>
#include <algorithm>

template <typename T>
std::vector<size_t> argsort(const std::vector<T> &v) {
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    std::stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

/*
  Selects elements of v using the passed in index.

  As an example, select(v, indices) is equivalent to 
  numpy's v[indices].
 */
template <typename T>
std::vector<T> select(const std::vector<T> &v, const std::vector<size_t> indices) {
    std::vector<T> w(indices.size());
    for (size_t i = 0; i < w.size(); i++) {
        w[i] = v[indices[i]];
    }
    return w;
}

#endif
