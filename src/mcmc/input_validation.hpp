// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "paratoric/types/types.hpp"

#include <stdexcept>

namespace paratoric::detail {

inline void validate_lattice_spec(const LatSpec& spec) {
    if (spec.basis != 'x' && spec.basis != 'z') {
        throw std::invalid_argument("basis must be 'x' or 'z'");
    }
    if (spec.lattice_type != "square"
        && spec.lattice_type != "cubic"
        && spec.lattice_type != "honeycomb"
        && spec.lattice_type != "triangular"
        && spec.lattice_type != "kagome") {
        throw std::invalid_argument("lattice_type is not supported");
    }
    if (spec.system_size <= 0) {
        throw std::invalid_argument("system_size must be a strictly positive integer");
    }
    if (!(spec.beta > 0.0)) {
        throw std::invalid_argument("beta must be strictly positive");
    }
    if (spec.boundaries != "periodic" && spec.boundaries != "open") {
        throw std::invalid_argument("boundaries must be 'periodic' or 'open'");
    }
    if (spec.default_spin != -1 && spec.default_spin != 1) {
        throw std::invalid_argument("default_spin must be -1 or 1");
    }
}

inline void validate_common_sim_spec(const SimSpec& spec) {
    if (spec.N_thermalization < 0) {
        throw std::invalid_argument("N_thermalization must be a non-negative integer");
    }
    if (spec.N_resamples <= 0) {
        throw std::invalid_argument("N_resamples must be a strictly positive integer");
    }
}

inline void validate_sampling_sim_spec(const SimSpec& spec) {
    validate_common_sim_spec(spec);
    if (spec.N_samples <= 0) {
        throw std::invalid_argument("N_samples must be a strictly positive integer");
    }
    if (spec.N_between_samples < 0) {
        throw std::invalid_argument("N_between_samples must be a non-negative integer");
    }
}

inline void validate_thermalization_config(const Config& config) {
    validate_lattice_spec(config.lat_spec);
    validate_common_sim_spec(config.sim_spec);
}

inline void validate_sample_config(const Config& config) {
    validate_lattice_spec(config.lat_spec);
    validate_sampling_sim_spec(config.sim_spec);
}

inline void validate_hysteresis_config(const Config& config) {
    validate_sample_config(config);

    if (config.param_spec.h_hys.empty()) {
        throw std::invalid_argument("h_hys must be non-empty");
    }
    if (config.param_spec.lmbda_hys.empty()) {
        throw std::invalid_argument("lmbda_hys must be non-empty");
    }
    if (config.param_spec.h_hys.size() != config.param_spec.lmbda_hys.size()) {
        throw std::invalid_argument("h_hys and lmbda_hys must have the same length");
    }
    if (config.out_spec.save_snapshots
        && config.out_spec.paths_out.size() != config.param_spec.h_hys.size()) {
        throw std::invalid_argument(
            "paths_out must match the hysteresis schedule length when save_snapshots is enabled"
        );
    }
}

} // namespace paratoric::detail
