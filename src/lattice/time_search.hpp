// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

/** @file
 * @brief Bound searches for sorted imaginary-time event histories.
 */
#pragma once

#include <iterator>

namespace paratoric::detail {

/**
 * @brief Find a lower or upper bound in a sorted event history.
 * @tparam Upper False finds the first value >= time; true finds the first > time.
 * @pre [first, last) is sorted in ascending order and contains no NaN values.
 * @return The bound iterator, or last when no value satisfies the bound.
 */

// Keep the number of remaining comparisons independent of the query value.
// The conditional iterator selection can use conditional-move instructions,
// avoiding an unpredictable branch at each level for random proposal times.
template<bool Upper, std::random_access_iterator Iterator>
inline Iterator time_bound(Iterator first, Iterator last, double time) {
    auto count = last - first;
    if (count == 0) return first;
    // A complete scan of a short contiguous history can compare multiple
    // timestamps per instruction, with independent loads instead of the
    // dependent loads of binary search. Retain logarithmic search for long
    // histories and avoid scan setup for tiny ones.
    if (count >= 8 && count <= 64) {
        decltype(count) rank = 0;
        for (decltype(count) i = 0; i < count; ++i)
            rank += Upper ? !(time < first[i]) : (first[i] < time);
        return first + rank;
    }
    while (count > 1) {
        const auto half = count / 2;
        const auto middle = first + half;
        const bool before = Upper ? !(time < *middle) : (*middle < time);
        first = before ? middle : first;
        count -= half;
    }
    const bool before = Upper ? !(time < *first) : (*first < time);
    return first + before;
}

/** @brief Return the first event >= time, or last. See time_bound(). */
template<std::random_access_iterator Iterator>
inline Iterator time_lower_bound(Iterator first, Iterator last, double time) {
    return time_bound<false>(first, last, time);
}

/** @brief Return the first event > time, or last. See time_bound(). */
template<std::random_access_iterator Iterator>
inline Iterator time_upper_bound(Iterator first, Iterator last, double time) {
    return time_bound<true>(first, last, time);
}

} // namespace paratoric::detail
