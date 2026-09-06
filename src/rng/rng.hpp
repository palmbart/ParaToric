// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <random>
#include <cstdint>

namespace paratoric::rng {

/**
 * @brief Uniform random bit generator wrapping std::mt19937_64.
 *
 * Default construction and copying reseed from std::random_device. Use set_seed()
 * for reproducible draws. Sharing an RNG via shared_ptr preserves one common
 * stream; it does not invoke the reseeding copy constructor.
 */
struct RNG {
    using result_type = std::mt19937_64::result_type;

    std::uint64_t seed_;
    std::mt19937_64 rng;

    explicit RNG(std::uint64_t s) : seed_(s), rng(seed_) {}

    RNG() : RNG(std::random_device{}()) {}

    result_type operator()() {
        return rng();
    }

    static constexpr result_type min() { return std::mt19937_64::min(); }
    static constexpr result_type max() { return std::mt19937_64::max(); }

    // Copying starts a fresh stream rather than duplicating the source state.
    RNG(RNG const&) : RNG() {}
    RNG& operator=(RNG const&) {
        seed_ = std::random_device{}();
        rng.seed(seed_);
        return *this;
    }

    /** @brief Restart the stream at seed s. */
    void set_seed(std::uint64_t s) {
        seed_ = s; rng.seed(seed_);
    }
    std::uint64_t get_seed() const { return seed_; }
};

/**
 * @brief Draw an unbiased integer in [0, bound).
 * @return Zero without consuming RNG state when bound is zero.
 * @note Multiply-and-reject uses 128-bit arithmetic where available; other
 *       platforms use std::uniform_int_distribution.
 */
inline std::uint64_t uniform_index(RNG& rng, std::uint64_t bound) {
    if (bound == 0) {
        return 0;
    }

#if defined(__SIZEOF_INT128__)
    while (true) {
        const auto product = static_cast<unsigned __int128>(rng()) * bound;
        const auto low = static_cast<std::uint64_t>(product);
        // threshold = 2^64 mod bound is strictly less than bound. For the
        // usual small lattice/event counts, low >= bound on almost every draw,
        // so the division is only needed in the rare rejection region.
        if (low >= bound || low >= -bound % bound) {
            return static_cast<std::uint64_t>(product >> 64);
        }
    }
#else
    std::uniform_int_distribution<std::uint64_t> dist(0, bound - 1);
    return dist(rng);
#endif
}

} // namespace paratoric::rng
