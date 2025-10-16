// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <random>
#include <cstdint>

namespace paratoric::rng {

// This RNG is an extension of std::mt19337_64 with additional features like manually setting seeds and a copy-costructor which sets a new random seed
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

    // copy = reseed
    RNG(RNG const&) : RNG() {}
    RNG& operator=(RNG const&) {
        seed_ = std::random_device{}();
        rng.seed(seed_);
        return *this;
    }

    void set_seed(std::uint64_t s) {
        seed_ = s; rng.seed(seed_);
    }
    std::uint64_t get_seed() const { return seed_; }
};

} // namespace paratoric::rng
