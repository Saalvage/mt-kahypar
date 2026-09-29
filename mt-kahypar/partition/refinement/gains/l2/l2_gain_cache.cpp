/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/gains/l2/l2_gain_cache.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_vector.h>

#include "mt-kahypar/datastructures/synchronized_edge_update.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/quadratic_delta.h"

namespace mt_kahypar {

template<typename PartitionedHypergraph>
void L2GainCache::initializeGainCache(const PartitionedHypergraph& partitioned_hg) {
  ASSERT(!_is_initialized, "Gain cache is already initialized");
  ASSERT(_k <= 0 || _k >= partitioned_hg.k(),
    "Gain cache was already initialized for a different k" << V(_k) << V(partitioned_hg.k()));
  allocateGainTable(partitioned_hg.topLevelNumNodes(), partitioned_hg.k());

  if constexpr (requires { partitioned_hg.partSumCutEdgeWeights(); }) {
    _edge_weight_sum_per_partition = &partitioned_hg.partSumCutEdgeWeights();
  }

  // TODO: Optimize.
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), partitioned_hg.initialNumNodes()),
    [&](tbb::blocked_range<HypernodeID>& r) {
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( partitioned_hg.nodeIsEnabled(u)) {
		  _gain_cache[penalty_index(u)].store(recomputePenaltyTerm(partitioned_hg, u), std::memory_order_relaxed);
          for (PartitionID to = 0; to < _k; ++to) {
            _gain_cache[benefit_index(u, to)].store(recomputeBenefitTerm(partitioned_hg, u, to), std::memory_order_relaxed);
          }
        }
      }
    });

  _is_initialized = true;
}

bool L2GainCache::triggersDeltaGainUpdate(const SynchronizedEdgeUpdate& sync_update) {
    return sync_update.pin_count_in_from_part_after == 0 ||
        sync_update.pin_count_in_from_part_after == 1 ||
        sync_update.pin_count_in_to_part_after == 1 ||
        sync_update.pin_count_in_to_part_after == 2 ||
        sync_update.pin_count_in_from_part_after == sync_update.edge_size - 1 ||
        sync_update.pin_count_in_from_part_after == sync_update.edge_size - 2 ||
        sync_update.pin_count_in_to_part_after == sync_update.edge_size ||
        sync_update.pin_count_in_to_part_after == sync_update.edge_size - 1;
}


template<typename PartitionedHypergraph>
void L2GainCache::deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                                   const SynchronizedEdgeUpdate& sync_update) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    const HypernodeID edge_size = sync_update.edge_size;
    if (edge_size == 1) {
        return;
    }

    const HyperedgeID he = sync_update.he;
    const PartitionID from = sync_update.from;
    const PartitionID to = sync_update.to;
    const HyperedgeWeight edge_weight = sync_update.edge_weight;
    const HypernodeID pin_count_in_from_part_after = sync_update.pin_count_in_from_part_after;
    const HypernodeID pin_count_in_to_part_after = sync_update.pin_count_in_to_part_after;

    if (pin_count_in_from_part_after == 1) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, from));
            if (partitioned_hg.partID(u) == from) {
                _gain_cache[penalty_index(u)].fetch_sub(edge_weight, std::memory_order_relaxed);
            }
        }
    } else if (pin_count_in_from_part_after == 0) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, from));
            _gain_cache[benefit_index(u, from)].fetch_add(edge_weight, std::memory_order_relaxed);
        }
        (*sync_update.edge_weight_sum_per_partition)[from].fetch_sub(edge_weight, std::memory_order_relaxed);
    }

    if (pin_count_in_to_part_after == 1) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, to));
            _gain_cache[benefit_index(u, to)].fetch_sub(edge_weight, std::memory_order_relaxed);
        }
        (*sync_update.edge_weight_sum_per_partition)[to].fetch_add(edge_weight, std::memory_order_relaxed);
    } else if (pin_count_in_to_part_after == 2) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, to));
            if (partitioned_hg.partID(u) == to) {
                _gain_cache[penalty_index(u)].fetch_add(edge_weight, std::memory_order_relaxed);
            }
        }
    }

    if (pin_count_in_from_part_after == edge_size - 1) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, from));
            _gain_cache[penalty_index(u)].fetch_sub(edge_weight, std::memory_order_relaxed);
            _gain_cache[benefit_index(u, from)].fetch_sub(edge_weight, std::memory_order_relaxed);
        }
        (*sync_update.edge_weight_sum_per_partition)[from].fetch_add(edge_weight, std::memory_order_relaxed);
    } else if (pin_count_in_from_part_after == edge_size - 2) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, from));
            _gain_cache[benefit_index(u, from)].fetch_add(edge_weight, std::memory_order_relaxed);
        }
    }

    if (pin_count_in_to_part_after == edge_size) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, to));
            _gain_cache[penalty_index(u)].fetch_add(edge_weight, std::memory_order_relaxed);
            _gain_cache[benefit_index(u, to)].fetch_add(edge_weight, std::memory_order_relaxed);
        }
        (*sync_update.edge_weight_sum_per_partition)[to].fetch_sub(edge_weight, std::memory_order_relaxed);
    } else if (pin_count_in_to_part_after == edge_size - 1) {
        for (const HypernodeID& u : partitioned_hg.pins(he)) {
            ASSERT(nodeGainAssertions(u, to));
            _gain_cache[benefit_index(u, to)].fetch_sub(edge_weight, std::memory_order_relaxed);
        }
    }
}

template<typename PartitionedHypergraph>
void L2GainCache::uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                                const HypernodeID u,
                                                const HypernodeID v,
                                                const HyperedgeID he,
                                                const HypernodeID pin_count_in_part_after) {
  throw 2;
}

template<typename PartitionedHypergraph>
void L2GainCache::uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                                    const HypernodeID u,
                                                    const HypernodeID v,
                                                    const HyperedgeID he) {
  throw 2;
}

template<typename PartitionedHypergraph>
void L2GainCache::initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                                   const HypernodeID u,
                                                   vec<Gain>& benefit_aggregator) {
  throw 2;
}


namespace {
#define L2_INITIALIZE_GAIN_CACHE(X) void L2GainCache::initializeGainCache(const X&)
#define L2_DELTA_GAIN_UPDATE(X) void L2GainCache::deltaGainUpdate(const X&,                     \
                                                                    const SynchronizedEdgeUpdate&)
#define L2_RESTORE_UPDATE(X) void L2GainCache::uncontractUpdateAfterRestore(const X&,          \
                                                                              const HypernodeID, \
                                                                              const HypernodeID, \
                                                                              const HyperedgeID, \
                                                                              const HypernodeID)
#define L2_REPLACEMENT_UPDATE(X) void L2GainCache::uncontractUpdateAfterReplacement(const X&,            \
                                                                                      const HypernodeID,   \
                                                                                      const HypernodeID,   \
                                                                                      const HyperedgeID)
#define L2_INIT_GAIN_CACHE_ENTRY(X) void L2GainCache::initializeGainCacheEntryForNode(const X&,           \
                                                                                        const HypernodeID,  \
                                                                                        vec<Gain>&)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(L2_INITIALIZE_GAIN_CACHE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(L2_DELTA_GAIN_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(L2_RESTORE_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(L2_REPLACEMENT_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(L2_INIT_GAIN_CACHE_ENTRY)

}  // namespace mt_kahypar
