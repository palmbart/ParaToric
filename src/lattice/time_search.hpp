// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

/** @file
 * @brief Geometry-independent imaginary-time histories and spin integrals.
 *
 * Histories are sorted ascending, contain finite times within one period, and
 * may contain coincident events (each occurrence toggles the spin). These
 * helpers do not own spins, select operator channels, maintain energy caches,
 * or synchronize related histories. Those responsibilities belong to the model.
 */
#pragma once

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <cstddef>
#include <format>
#include <functional>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

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

/**
 * @brief Borrow a spin and its history, independently of site or bond storage.
 * @pre spin is +1 or -1.
 * @note spin is the incoming value before any events at zero. Reacquire this
 *       view after mutating the spin or its history. Normally supply all flips;
 *       reduced histories require a model-specific proof of equivalence.
 */
struct WorldlineView {
    int spin;
    std::span<const double> flips;
};

/** @brief First exact event rank; throw std::runtime_error if absent. */
inline std::size_t time_index(std::span<const double> times, double tau) {
    const auto it = time_lower_bound(times.begin(), times.end(), tau);
    if (it == times.end() || *it != tau) [[unlikely]] {
        throw std::runtime_error(std::format("time_index: There is no event at {}.", tau));
    }
    return static_cast<std::size_t>(it - times.begin());
}

/** @brief Next event strictly after tau, wrapping; return tau for an empty history. */
inline double next_event_time(std::span<const double> times, double tau) {
    if (times.empty()) [[unlikely]] return tau;
    const auto it = time_upper_bound(times.begin(), times.end(), tau);
    return it != times.end() ? *it : times.front();
}

/** @brief Previous event strictly before tau, wrapping; return tau for an empty history. */
inline double previous_event_time(std::span<const double> times, double tau) {
    if (times.empty()) [[unlikely]] return tau;
    const auto it = time_lower_bound(times.begin(), times.end(), tau);
    return it != times.begin() ? *(it - 1) : times.back();
}

/** @brief Insert after any equal-time events. Invalidates borrowed history views. */
inline void insert_time(std::vector<double>& times, double tau) {
    if (times.empty()) [[unlikely]] {
        times.emplace_back(tau);
    } else [[likely]] {
        times.insert(time_upper_bound(times.begin(), times.end(), tau), tau);
    }
}

/** @brief Insert ordered events (left <= right), preserving coincident events. */
inline void insert_time_pair(std::vector<double>& times, double left, double right) {
    if (right < left) [[unlikely]] {
        throw std::invalid_argument("insert_time_pair: event times must be ordered.");
    }
    if (times.empty()) [[unlikely]] {
        times.push_back(left);
        times.push_back(right);
    } else [[likely]] {
        auto it = time_upper_bound(times.begin(), times.end(), right);
        it = times.insert(it, right);
        times.insert(time_lower_bound(times.begin(), it, left), left);
    }
}

/** @brief Erase the first exact occurrence; throw std::runtime_error if absent. */
inline void erase_time(std::vector<double>& times, double tau) {
    times.erase(times.begin() + time_index(times, tau));
}

/**
 * @brief Erase ordered events (left < right).
 * @note Erases the first right occurrence, then the last left occurrence before
 *       it. If left is absent, right has already been erased. Callers must
 *       ensure both events exist in every history they intend to synchronize.
 */
inline void erase_time_pair(std::vector<double>& times, double left, double right) {
    if (right < left) [[unlikely]] {
        throw std::invalid_argument("erase_time_pair: event times must be ordered.");
    }
    const auto it = times.erase(times.begin() + time_index(times, right));
    const auto reverse = std::find(std::make_reverse_iterator(it), times.rend(), left);
    if (reverse == times.rend()) [[unlikely]] {
        throw std::runtime_error(std::format("erase_time_pair: There is no event at {}.", left));
    }
    times.erase(std::prev(reverse.base()));
}

/**
 * @brief Move an indexed event within its neighboring-event window.
 * @pre index is valid; no other event is crossed except through the time origin.
 * @param no_wrap False moves the first event to the end or the last to the front.
 * @note The caller must also reverse the origin spin when crossing the origin.
 */
inline void move_time(std::vector<double>& times, std::size_t index, double tau, bool no_wrap) {
    if (no_wrap) [[likely]] {
        times[index] = tau;
    } else if (index == 0) {
        times.front() = tau;
        std::rotate(times.begin(), times.begin() + 1, times.end());
    } else {
        std::rotate(times.rbegin(), times.rbegin() + 1, times.rend());
        times.front() = tau;
    }
}

/**
 * @brief Shift the time origin to tau, retaining sorted event order.
 * @pre beta > 0, 0 <= tau < beta, event times in [0, beta).
 * @param scratch Reusable buffer, distinct from times, grown only as needed.
 * @return Number of events strictly before tau. Reverse the origin spin if odd.
 * @note Events exactly at the cut become epsilon, preserving the existing ETC
 *       convention that they occur just after the new origin.
 */
inline std::size_t rotate_times(
    std::vector<double>& times, double tau, double beta, std::vector<double>& scratch
) {
    const auto pivot = static_cast<std::size_t>(
        time_lower_bound(times.begin(), times.end(), tau) - times.begin());
    if (scratch.size() < times.size()) scratch.resize(times.size());
    const auto shift = [&](double t) {
        t -= tau;
        if (t < 0.) t += beta;
        else if (t >= beta) t -= beta;
        if (t == 0.) t += std::numeric_limits<double>::epsilon();
        return t;
    };
    for (std::size_t i = pivot; i < times.size(); ++i) {
        scratch[i - pivot] = shift(times[i]);
    }
    for (std::size_t i = 0; i < pivot; ++i) {
        scratch[times.size() - pivot + i] = shift(times[i]);
    }
    std::copy_n(scratch.begin(), times.size(), times.begin());
    return pivot;
}

/** @brief Spin after all events at tau; does not wrap the query at beta. */
inline int spin_at_time(WorldlineView line, double tau) {
    const auto hi = time_upper_bound(line.flips.begin(), line.flips.end(), tau);
    return ((hi - line.flips.begin()) & 1) ? -line.spin : line.spin;
}

/**
 * @brief Integrate a spin times a weight over an ordered, non-wrapping interval.
 * @param weight_integral Callable returning the weight integral on [a, b].
 * @return Bare integral, without coupling or Hamiltonian sign; zero if left >= right.
 */
template<class WeightIntegral>
inline double weighted_spin_integral(
    WorldlineView line, double left, double right, WeightIntegral weight_integral
) {
    if (left >= right) return 0.0;
    const auto& flips = line.flips;
    const auto lo = time_lower_bound(flips.begin(), flips.end(), left);
    int spin = ((lo - flips.begin()) & 1) ? -line.spin : line.spin;
    double integral = 0.0;
    double prev = left;
    for (auto it = lo; it != flips.end() && *it <= right; ++it) {
        const double curr = *it;
        integral += spin * weight_integral(prev, curr);
        spin = -spin;
        prev = curr;
    }
    if (prev < right) integral += spin * weight_integral(prev, right);
    return integral;
}

/** @brief Bare spin integral; ordered, non-wrapping bounds, zero if left >= right. */
inline double spin_integral(WorldlineView line, double left, double right) {
    return weighted_spin_integral(line, left, right, [](double a, double b) { return b - a; });
}

/**
 * @brief Bare spin integral on an interval with no events strictly inside it.
 * @pre left < right. A supplied rank/time pair must identify the first occurrence
 *       of known_time in the history.
 */
inline double spin_integral_no_inner_flips(
    WorldlineView line, double left, double right, int known_index = -1, double known_time = 0.0
) {
    if (known_index >= 0 && right == known_time) {
        const int spin = (known_index & 1) ? -line.spin : line.spin;
        return (right - left) * spin;
    }
    auto it = known_index >= 0 && left == known_time
        ? line.flips.begin() + known_index
        : time_lower_bound(line.flips.begin(), line.flips.end(), left);
    int spin = ((it - line.flips.begin()) & 1) ? -line.spin : line.spin;
    while (it != line.flips.end() && *it == left) {
        spin = -spin;
        ++it;
    }
    return (right - left) * spin;
}

/** @brief Antiderivative of min(tau, beta - tau), clamped outside [0, beta]. */
inline double triangular_weight_primitive(double tau, double beta) {
    if (tau <= 0.0) return 0.0;
    if (tau <= 0.5 * beta) return 0.5 * tau * tau;
    if (tau <= beta) return beta * tau - 0.5 * tau * tau - 0.25 * beta * beta;
    return 0.25 * beta * beta;
}

/**
 * @brief Integrate a spin product by merging histories with equal-time parity.
 * @param lines Worldline views, or model-owned entries with a get_worldline adapter.
 * @param get_worldline Maps each entry to WorldlineView; does not copy histories.
 * @pre Ordered, non-wrapping bounds within the stored period.
 * @return Bare integral; zero for reversed/empty intervals. An empty product is 1.
 * @note Uses the supplied histories in full; choosing a reduced operator history
 *       is the model's responsibility. Small products retain inline merge storage.
 */
template<std::ranges::sized_range Range, class GetWorldline = std::identity>
inline double spin_product_integral(
    const Range& lines, double imag_time_1, double imag_time_2, GetWorldline get_worldline = {}
) {
    if (imag_time_1 >= imag_time_2) {
        return 0.0;
    }

    // Merge histories backwards, combining equal-time flips by parity. This
    // avoids allocating and sorting a combined list for every tuple integral.
    struct ReverseStream {
        std::span<const double>::iterator begin;
        std::span<const double>::iterator it; // one past current event
    };

    boost::container::small_vector<ReverseStream, 8> streams;
    streams.reserve(std::ranges::size(lines));

    int spin_at_t2 = 1;
    for (const auto& entry : lines) {
        const WorldlineView line = get_worldline(entry);
        const auto spin_flips = line.flips;
        auto hi = time_upper_bound(spin_flips.begin(), spin_flips.end(), imag_time_2);

        int s = line.spin;
        if (((hi - spin_flips.begin()) & 1) != 0) {
            s = -s;
        }
        spin_at_t2 *= s;

        if (hi != spin_flips.begin() && *(hi - 1) >= imag_time_1) {
            streams.push_back({spin_flips.begin(), hi});
        }
    }

    double odd_len = 0.0;
    bool odd = false;
    double t_next = imag_time_2;

    while (!streams.empty()) {
        double max_t = *(streams[0].it - 1);
        for (size_t i = 1; i < streams.size(); ++i) {
            const double cand = *(streams[i].it - 1);
            if (cand > max_t) {
                max_t = cand;
            }
        }
        if (odd) {
            odd_len += t_next - max_t;
        }

        bool toggles = false;
        size_t i = 0;
        while (i < streams.size()) {
            auto& stream = streams[i];
            if (*(stream.it - 1) != max_t) {
                ++i;
                continue;
            }

            toggles = !toggles;
            --stream.it;

            if (stream.it == stream.begin || *(stream.it - 1) < imag_time_1) {
                streams[i] = streams.back();
                streams.pop_back();
                continue;
            }
            ++i;
        }

        if (toggles) {
            odd = !odd;
        }
        t_next = max_t;
    }

    if (odd) {
        odd_len += t_next - imag_time_1;
    }

    const double interval_len = imag_time_2 - imag_time_1;
    return static_cast<double>(spin_at_t2) * (interval_len - 2.0 * odd_len);
}

} // namespace paratoric::detail
