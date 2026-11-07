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

/**
 * @brief Draw a Poisson count conditioned on its parity; mean must be positive.
 *
 * At small mean, odd counts use the size-biased proposal 1 + Poisson(mean),
 * accepting odd k with probability 1/k. Its accepted mass is proportional to
 * mean^k/k!, without a rejection cost proportional to 1/mean. For zero mean
 * (e.g. product underflow), return the limiting conditional count.
 */
inline std::uint64_t poisson_with_parity(RNG& rng, double mean, bool odd) {
    if (mean == 0.) return odd ? 1 : 0;
    std::poisson_distribution<std::uint64_t> poisson(mean);
    std::uniform_real_distribution<double> uniform(0., 1.);
    if (odd && mean < 1.) {
        while (true) {
            const auto k = 1 + poisson(rng);
            if (k % 2 == 1 && uniform(rng) < 1. / static_cast<double>(k)) return k;
        }
    }
    while (true) {
        const auto k = poisson(rng);
        if ((k % 2 == 1) == odd) return k;
    }
}

} // namespace paratoric::rng
