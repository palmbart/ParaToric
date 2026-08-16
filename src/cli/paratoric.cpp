// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#include "io/io.hpp"
#include "mcmc/input_validation.hpp"
#include "paratoric/types/types.hpp"

#include <boost/log/core.hpp> 
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp> 
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    try {
        paratoric::Config config{};
        std::string simulation;
        int process_index = 0;

        po::options_description desc("Run extended toric code QMC simulation");
        desc.add_options()
        ("help,h", "Display this help message.")
        ("simulation,sim", po::value(&simulation), "Simulation you want to run.")
        ("N_samples,Ns", po::value(&config.sim_spec.N_samples), "Number of samples taken during the simulation.")
        ("N_thermalization,Nth", po::value(&config.sim_spec.N_thermalization), "Number of thermalization steps.")
        ("N_between_samples,Nbs", po::value(&config.sim_spec.N_between_samples), "Number of steps between samples.")
        ("beta,bet", po::value(&config.lat_spec.beta), "Inverse temperature.")
        ("h_constant,hc", po::value(&config.param_spec.h), "Value of the constant h (electric field term).")
        ("h_hysteresis,hhys", po::value(&config.param_spec.h_hys)->multitoken(), "Hysteresis values of the constant h (electric field term).")
        ("h_constant_therm,hct", po::value(&config.param_spec.h_therm), "Value of the constant h for thermalization (electric field term).")
        ("mu_constant,muc", po::value(&config.param_spec.mu), "Value of the constant mu_constant (star term).")
        ("J_constant,Jc", po::value(&config.param_spec.J), "Value of the constant J (plaquette term).")
        ("lmbda_constant,lmbdac", po::value(&config.param_spec.lmbda), "Value of the constant lmbda (gauge field term).")
        ("lmbda_hysteresis,lmbdahys", po::value(&config.param_spec.lmbda_hys)->multitoken(), "Hysteresis values of the constant lmbda (gauge field term).")
        ("lmbda_constant_therm,lmbdact", po::value(&config.param_spec.lmbda_therm), "Value of the constant lmbda for thermalization (gauge field term).")
        ("N_resamples,Nr", po::value(&config.sim_spec.N_resamples), "Number of bootstrap resamples.")
        ("custom_therm,cth", po::value(&config.sim_spec.custom_therm), "Whether custom thermalization is used (to probe hysteresis).")
        ("observables,obs", po::value(&config.sim_spec.observables)->multitoken(), "Observables that are measured.") // can take multiple observables at once
        ("seed,s", po::value(&config.sim_spec.seed), "Seed for the pseudorandom number generator.")
        ("basis,bas", po::value(&config.lat_spec.basis), "Spin basis (\"x\" or \"z\").")
        ("output_directory,outdir", po::value(&config.out_spec.path_out), "Directory where the output is stored.")
        ("folder_name,fn", po::value(&config.out_spec.folder_name), "Name of the output subfolder.")
        ("folder_names,fns", po::value(&config.out_spec.folder_names)->multitoken(), "Directories where hysteresis output is stored.")
        ("snapshots,snap", po::value(&config.out_spec.save_snapshots), "Whether snapshots should be saved.")
        ("full_time_series,fts", po::value(&config.out_spec.full_time_series), "Whether full time series should be saved.")
        ("process_index,procid", po::value(&process_index), "Identifier of process (for debugging).")
        ("lattice_type,lat", po::value(&config.lat_spec.lattice_type), "Type of lattice used.")
        ("system_size,L", po::value(&config.lat_spec.system_size), "System size of lattice (one coordinate).")
        ("boundaries,bound", po::value(&config.lat_spec.boundaries), "Boundary condition of the lattice (periodic or open).")
        ("default_spin,dsp", po::value(&config.lat_spec.default_spin), "Default spin (electric field) for lattice initialization.");

        // With long-option disguise enabled, Boost otherwise treats "-h" as
        // an ambiguous prefix of the h-related options before resolving the
        // explicit short alias.
        for (int i = 1; i < argc; ++i) {
            if (std::string_view(argv[i]) == "-h") {
                std::cout << desc << '\n';
                return 0;
            }
        }

        po::variables_map vm;
        po::command_line_parser parser(argc, argv);
        parser.options(desc).style(
            po::command_line_style::default_style
            | po::command_line_style::allow_long_disguise
        );
        po::store(parser.run(), vm);
        po::notify(vm);

        if (vm.count("help") || argc == 1) {
            std::cout << desc << '\n';
            return 0;
        }
        if (!vm.count("simulation")) {
            throw std::invalid_argument("simulation must be specified");
        }

        if (simulation == "etc_sample") {
            paratoric::detail::validate_sample_config(config);
        } else if (simulation == "etc_hysteresis") {
            if (config.out_spec.folder_names.size() != config.param_spec.h_hys.size()) {
                throw std::invalid_argument(
                    "folder_names must match the hysteresis schedule length"
                );
            }
            auto config_to_validate = config;
            config_to_validate.out_spec.paths_out.resize(
                config.out_spec.folder_names.size()
            );
            paratoric::detail::validate_hysteresis_config(config_to_validate);
        } else if (simulation == "etc_thermalization") {
            paratoric::detail::validate_thermalization_config(config);
        } else {
            throw std::invalid_argument(
                "simulation must be 'etc_sample', 'etc_hysteresis', or 'etc_thermalization'"
            );
        }

        BOOST_LOG_TRIVIAL(info) << "Subprocess [ID=" << process_index 
                                << "] SPAWNED with [β=" << config.lat_spec.beta
                                << "], [μ=" << config.param_spec.mu
                                << "], [J=" << config.param_spec.J
                                << "], [h=" << config.param_spec.h
                                << "], [λ=" << config.param_spec.lmbda << "].";

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        if (simulation == "etc_sample") {
            auto io = paratoric::IO();
            io.etc_sample(config);
        } else if (simulation == "etc_hysteresis") {
            auto io = paratoric::IO();
            io.etc_hysteresis(config);
        } else if (simulation == "etc_thermalization") {
            auto io = paratoric::IO();
            io.etc_thermalization(config);
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        const auto dt = end - begin;
        const auto days = std::chrono::duration_cast<std::chrono::days>(dt);
        const std::chrono::hh_mm_ss<std::chrono::microseconds> hms{floor<std::chrono::microseconds>(dt - days)};

        BOOST_LOG_TRIVIAL(info)
        << std::format(
            "Subprocess [ID={}] COMPLETED in [{}-{:02}:{:02}:{:02}.{:06} h].",
            process_index,
            days.count(),
            hms.hours().count(),
            hms.minutes().count(),
            hms.seconds().count(),
            hms.subseconds().count());

        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "Exception caught: " << exc.what() << std::endl;
        return 2;
    } 
}
