#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {

struct delta_vector {
    template <typename T>
    delta_vector(const T* base) {
        _delta.clear();
        _delta.resize(base->size());
        for (std::size_t i = 0; i < base->size(); ++i) {
            if constexpr (requires{ (*base)[i].load(std::memory_order_relaxed); }) {
                _delta[i] = (*base)[i].load(std::memory_order_relaxed);
            } else {
                _delta[i] = (*base)[i];
            }
        }
    }

    template <typename PartitionedHypergraph>
    delta_vector(const PartitionedHypergraph& phg) {
        _delta.clear();
        _delta.resize(phg.k());
        for (std::size_t i = 0; i < _delta.size(); ++i) {
            _delta[i] = phg.partSumCutEdgeWeight(i);
        }
    }

    void add_at(std::size_t idx, HyperedgeWeight delta) {
        _delta[idx] += delta;
    }

    HyperedgeWeight max() {
        HyperedgeWeight max_delta = 0;
        for (HyperedgeWeight i : _delta) {
            max_delta = std::max(max_delta, i);
        }
        return max_delta;
    }

private:
    std::vector<HyperedgeWeight> _delta;
};

}
