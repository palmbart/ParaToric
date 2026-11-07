// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

namespace paratoric::detail {

/**
 * @brief One-event parity proposal for a Poisson expansion with mean beta*c.
 * From count k, choose birth with probability x/(x+k), or uniformly delete an
 * existing event. The Hastings factor, after cancellation of expansion powers,
 * is (x+k)/(x+k'), where k' = k +/- 1. Logarithms retain tiny-coupling limits.
 */
struct PoissonParityProposal {
    double mean;
    double zero_to_one_log;

    PoissonParityProposal(double beta, double coupling)
        : mean(beta * coupling),
          zero_to_one_log((mean == 0. ? std::log(beta) + std::log(coupling) : std::log(mean))
                          - std::log1p(mean)) {}

    double birth_probability(std::size_t count) const {
        return count == 0 ? 1. : mean / (mean + static_cast<double>(count));
    }

    double log_ratio(std::size_t count, bool create) const {
        if (count == 0) return zero_to_one_log;
        if (count == 1 && !create) return -zero_to_one_log;
        const auto proposed_count = create ? count + 1 : count - 1;
        return std::log1p((create ? -1. : 1.) / (mean + static_cast<double>(proposed_count)));
    }
};

/**
 * @brief Kernel basis over GF(2) of the edge/tuple incidence matrix.
 * @param columns Edge indices flipped by each tuple, in tuple-index order.
 * @param edge_count Number of rows of the incidence matrix.
 * @return Independent lists of tuple indices whose product is the identity.
 *
 * Sparse column elimination tracks the original columns in each pivot. Every
 * dependent column supplies one relation with a distinct last tuple index, so
 * the returned relations span all trace-compatible tuple-count parities.
 */
inline std::vector<std::vector<int>> tuple_parity_kernel(
    const std::vector<std::vector<int>>& columns, int edge_count
) {
    std::vector<std::vector<int>> pivots(edge_count), origins(edge_count), kernel;
    const auto xor_into = [](std::vector<int>& a, const std::vector<int>& b) {
        std::vector<int> result;
        result.reserve(a.size() + b.size());
        std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(),
                                    std::back_inserter(result));
        a = std::move(result);
    };

    for (std::size_t i = 0; i < columns.size(); ++i) {
        auto support = columns[i];
        std::sort(support.begin(), support.end());
        // Repeated occurrences of an edge cancel rather than becoming one flip.
        std::vector<int> odd_support;
        for (std::size_t j = 0; j < support.size();) {
            auto next = j + 1;
            while (next < support.size() && support[next] == support[j]) ++next;
            if ((next - j) % 2 != 0) odd_support.push_back(support[j]);
            j = next;
        }
        support = std::move(odd_support);
        std::vector<int> combination{static_cast<int>(i)};
        while (!support.empty() && !pivots[support.front()].empty()) {
            const int pivot = support.front();
            xor_into(support, pivots[pivot]);
            xor_into(combination, origins[pivot]);
        }
        if (support.empty()) {
            kernel.push_back(std::move(combination));
        } else {
            const int pivot = support.front();
            pivots[pivot] = std::move(support);
            origins[pivot] = std::move(combination);
        }
    }
    return kernel;
}

/**
 * @brief Replace a kernel basis by independent short relations, retaining its span.
 * Candidates must already be identity products. Keep their original short support
 * and supplement them with the remaining kernel directions. Reduced supplements
 * can be shorter wrapping surfaces than the original elimination relations.
 */
inline std::vector<std::vector<int>> prefer_local_tuple_relations(
    const std::vector<std::vector<int>>& kernel,
    const std::vector<std::vector<int>>& local, int tuple_count
) {
    std::vector<std::vector<int>> pivots(tuple_count), result;
    result.reserve(kernel.size());
    const auto add = [&](const std::vector<int>& relation, bool keep_original) {
        auto reduced = relation;
        std::sort(reduced.begin(), reduced.end());
        while (!reduced.empty() && !pivots[reduced.front()].empty()) {
            std::vector<int> next;
            const auto& pivot = pivots[reduced.front()];
            next.reserve(reduced.size() + pivot.size());
            std::set_symmetric_difference(reduced.begin(), reduced.end(), pivot.begin(), pivot.end(),
                                         std::back_inserter(next));
            reduced = std::move(next);
        }
        if (!reduced.empty()) {
            result.push_back(keep_original || relation.size() < reduced.size() ? relation : reduced);
            pivots[reduced.front()] = std::move(reduced);
        }
    };
    for (const auto& relation : local) add(relation, true);
    for (const auto& relation : kernel) add(relation, false);
    return result;
}

} // namespace paratoric::detail
