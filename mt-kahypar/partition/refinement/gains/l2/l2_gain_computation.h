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

#include <vector>

#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/partition/refinement/gains/l2/l2_attributed_gains.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

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
      PartitionID connectivity = phg.connectivity(he);
      HypernodeID pin_count_in_from_part = phg.pinCountInPart(he, from);
      HyperedgeWeight weight = phg.edgeWeight(he);
      HypernodeID edge_size = phg.edgeSize(he);
      if (connectivity == 1 && edge_size > 1) {
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        HyperedgeWeight part_sum_weight = phg.partSumCutEdgeWeight(from);
        HyperedgeWeight part_sum_weight_after = part_sum_weight + weight;
        isolated_block_gain += part_sum_weight_after * part_sum_weight_after - part_sum_weight * part_sum_weight;
      } else if (pin_count_in_from_part == 1) {
        HyperedgeWeight part_sum_weight = phg.partSumCutEdgeWeight(from);
        HyperedgeWeight part_sum_weight_after = part_sum_weight - weight;
        HyperedgeWeight score_gain = part_sum_weight_after * part_sum_weight_after - part_sum_weight * part_sum_weight;
        for (const PartitionID& to : phg.connectivitySet(he)) {
          // In case there are only two blocks contained in the current
          // hyperedge and only one pin left in the from part of the hyperedge,
          // we would make the current hyperedge a non-cut hyperedge when moving
          // vertex hn to the other block.
          if (from != to) {
            tmp_scores[to] -= score_gain;

            // Moving to this block would make it a non-cut hyperedge.
            if (phg.pinCountInPart(he, to) == edge_size - 1) {
              HyperedgeWeight part_sum_weight_to = phg.partSumCutEdgeWeight(to);
              HyperedgeWeight part_sum_weight_to_after = part_sum_weight_to - weight;
              tmp_scores[to] -= part_sum_weight_to_after * part_sum_weight_to_after - part_sum_weight_to * part_sum_weight_to;
            }
          }
        }
      }

      // TODO: Possibly optimize.
      for (const PartitionID& to : phg.connectivitySet(he)) {
        if (from != to && phg.pinCountInPart(he, to) == 0) {
          HyperedgeWeight part_sum_weight_to = phg.partSumCutEdgeWeight(to);
          HyperedgeWeight part_sum_weight_to_after = part_sum_weight_to + weight;
          tmp_scores[to] -= part_sum_weight_to_after * part_sum_weight_to_after - part_sum_weight_to * part_sum_weight_to;
        }
      }
    }
  }

  // ! Computes only the gain (cut increase) for moving out of the current block.
  template<typename PartitionedHypergraph>
  static Gain computeIsolatedBlockGain(const PartitionedHypergraph& phg, const HypernodeID hn) {
    Gain isolated_block_sum_cut_edge_weight = 0;
    Gain existing_block_initial_sum_cut_edge_weight = phg.partSumCutEdgeWeight(phg.partID(hn));
    Gain existing_block_sum_cut_edge_weight = existing_block_initial_sum_cut_edge_weight;
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
      if (phg.edgeSize(he) == 1) continue;

      auto weight = phg.edgeWeight(he);
      isolated_block_sum_cut_edge_weight += weight;

      PartitionID connectivity = phg.connectivity(he);
      if (connectivity == 1) {
        // In case, the hyperedge is a non-cut hyperedge, we would increase
        // the cut, if we move vertex hn to an other block.
        existing_block_sum_cut_edge_weight += phg.edgeWeight(he);
      }
    }
    Gain existing_block_contribution = existing_block_sum_cut_edge_weight * existing_block_sum_cut_edge_weight
      - existing_block_initial_sum_cut_edge_weight * existing_block_initial_sum_cut_edge_weight;
    return existing_block_contribution + isolated_block_sum_cut_edge_weight * isolated_block_sum_cut_edge_weight;
  }

  HyperedgeWeight gain(const Gain to_score,
                       const Gain isolated_block_gain) {
    return -to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID) {
    // Do nothing
  }
};
}  // namespace mt_kahypar
