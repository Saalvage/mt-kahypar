/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#pragma once

#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/partition/refinement/gains/l2/l2_attributed_gains.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/utils/quadratic_delta.h"

namespace mt_kahypar {

class L2GainComputation : public GainComputationBase<L2GainComputation, L2AttributedGains> {
  using Base = GainComputationBase<L2GainComputation, L2AttributedGains>;

  static constexpr bool enable_heavy_assert = false;

 public:
  using RatingMap = typename Base::RatingMap;

  static constexpr bool is_independent_of_block = true;

  L2GainComputation(const Context& context,
                     bool disable_randomization = false) :
    Base(context, disable_randomization) { }

  // The concept of first computing an isolated block gain does not work for the L2 metric, as the isolated block would
  // contribute a new addend to the total metric which is never equal to a move to a non-empty block.
  template<typename PartitionedHypergraph>
  void precomputeGains(const PartitionedHypergraph& phg,
                       const HypernodeID hn,
                       RatingMap& tmp_scores,
                       Gain& isolated_block_gain,
                       const bool) {
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");
    PartitionID from = phg.partID(hn);
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
      HypernodeID edge_size = phg.edgeSize(he);

      // Node is only part of net, can never cut.
      if (edge_size == 1) continue;

      HypernodeID pin_count_in_from_part = phg.pinCountInPart(he, from);
      HyperedgeWeight weight = phg.edgeWeight(he);
      // TODO: Use connectivity here?
      if (pin_count_in_from_part == edge_size) {
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        isolated_block_gain += weight;
      } else if (pin_count_in_from_part == 1) {
        // The current block would no longer be affected by the cut.
        isolated_block_gain -= weight;
      }

      for (const PartitionID& to : phg.connectivitySet(he)) {
        if (from == to) continue;

        HypernodeID pin_count_in_to_part = phg.pinCountInPart(he, to);
        if (pin_count_in_to_part == 0) {
          tmp_scores[to] -= weight;
        } else if (pin_count_in_to_part == edge_size - 1) {
          tmp_scores[to] += weight;
        }
      }
    }

    isolated_block_gain = quadratic_delta(phg.partSumCutEdgeWeight(from), isolated_block_gain);
    for (auto& [to, weight_delta] : tmp_scores) {
      weight_delta = quadratic_delta(phg.partSumCutEdgeWeight(to), weight_delta);
    }
  }

  // ! Computes only the gain (cut increase) for moving out of the current block.
  template<typename PartitionedHypergraph>
  static Gain computeIsolatedBlockGain(const PartitionedHypergraph& phg, const HypernodeID hn) {
    Gain isolated_block_sum_cut_edge_weight = 0;
    Gain existing_block_sum_cut_edge_weight = 0;
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
      if (phg.edgeSize(he) == 1) continue;

      auto weight = phg.edgeWeight(he);
      isolated_block_sum_cut_edge_weight += weight;

      if (phg.connectivity(he) == 1) {
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        existing_block_sum_cut_edge_weight += weight;
      } else if (phg.pinCountInPart(he, phg.partID(hn)) == 1) {
        // The current block would no longer be affected by the cut.
        existing_block_sum_cut_edge_weight -= weight;
      }
    }
    Gain existing_block_contribution = quadratic_delta(phg.partSumCutEdgeWeight(phg.partID(hn)), existing_block_sum_cut_edge_weight);
    return existing_block_contribution + isolated_block_sum_cut_edge_weight * isolated_block_sum_cut_edge_weight;
  }

  HyperedgeWeight gain(const Gain to_score,
                       const Gain isolated_block_gain) {
    return isolated_block_gain - to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID) {
    // Do nothing
  }
};
}  // namespace mt_kahypar
