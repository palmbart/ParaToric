// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#define BOOST_TEST_MODULE TestInputValidation

#include "paratoric/mcmc/extended_toric_code.hpp"
#include "paratoric/mcmc/extended_toric_code_c.h"
#include "paratoric/types/types.hpp"

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace {

paratoric::Config sample_config() {
    return paratoric::Config{
        .sim_spec = {
            .N_samples = 2,
            .N_thermalization = 0,
            .N_between_samples = 0,
            .N_resamples = 1,
            .observables = {}
        },
        .lat_spec = {
            .basis = 'x',
            .lattice_type = "square",
            .system_size = 2,
            .beta = 1.0,
            .boundaries = "periodic",
            .default_spin = 1
        }
    };
}

ptc_config_t c_sample_config() {
    ptc_config_t config{};
    config.sim.N_samples = 2;
    config.sim.N_thermalization = 0;
    config.sim.N_between_samples = 0;
    config.sim.N_resamples = 1;
    config.lat.basis = 'x';
    config.lat.lattice_type = "square";
    config.lat.system_size = 2;
    config.lat.beta = 1.0;
    config.lat.boundaries = "periodic";
    config.lat.default_spin = 1;
    return config;
}

} // namespace

BOOST_AUTO_TEST_CASE(sample_counts_are_validated) {
    auto config = sample_config();

    config.sim_spec.N_samples = 0;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.sim_spec.N_thermalization = -1;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.sim_spec.N_between_samples = -1;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.sim_spec.N_resamples = 0;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(lattice_options_are_validated) {
    auto config = sample_config();
    config.lat_spec.lattice_type = "unknown";
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.lat_spec.system_size = 0;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.lat_spec.beta = 0.0;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.lat_spec.boundaries = "unknown";
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.lat_spec.default_spin = 0;
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );

    config = sample_config();
    config.lat_spec.basis = 'y';
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_sample(config),
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(hysteresis_schedules_are_validated) {
    auto config = sample_config();
    config.param_spec.h_hys = {0.0, 0.1};
    config.param_spec.lmbda_hys = {0.0};

    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_hysteresis(config),
        std::invalid_argument
    );

    config.param_spec.lmbda_hys = {0.0, 0.1};
    config.out_spec.save_snapshots = true;
    config.out_spec.paths_out = {"first"};
    BOOST_CHECK_THROW(
        paratoric::ExtendedToricCode::get_hysteresis(config),
        std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(c_api_reports_invalid_values) {
    ptc_handle_t* handle = nullptr;
    BOOST_REQUIRE_EQUAL(ptc_create(&handle), PTC_STATUS_OK);

    auto config = c_sample_config();
    config.sim.N_samples = 0;
    ptc_result_t result{};

    BOOST_CHECK_EQUAL(
        ptc_get_sample(handle, &config, &result),
        PTC_STATUS_INVALID_ARGUMENT
    );
    BOOST_CHECK_NE(std::string(ptc_last_error()).find("N_samples"), std::string::npos);

    config = c_sample_config();
    config.sim.N_observables = 1;
    config.sim.observables = nullptr;
    BOOST_CHECK_EQUAL(
        ptc_get_sample(handle, &config, &result),
        PTC_STATUS_INVALID_ARGUMENT
    );
    BOOST_CHECK_NE(std::string(ptc_last_error()).find("observables"), std::string::npos);

    ptc_result_destroy(&result);
    ptc_destroy(handle);
}
