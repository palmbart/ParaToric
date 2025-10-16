// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "io/io.hpp"
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
#include <vector>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    try {
        std::string simulation;
        int N_samples;
        int N_thermalization;
        int N_between_samples;
        double beta;
        double h_constant;
        double h_constant_therm;
        std::vector<double> h_hys;
        double mu_constant;
        double J_constant;
        double lmbda_constant;
        double lmbda_constant_therm;
        std::vector<double> lmbda_hys;
        int N_resamples;
        bool custom_therm;
        std::vector<std::string> obs;
        int seed;
        char basis;
        std::string lattice_type;
        int system_size;
        std::string boundaries;
        int default_spin;
        std::filesystem::path path_output_directory;
        std::string folder_name;
        std::vector<std::string> folder_names;
        bool save_snapshots;
        bool full_time_series;
        int process_index;

        po::options_description desc("Run extended toric code QMC simulation");
        desc.add_options()
        ("simulation,sim", po::value(&simulation), "Simulation you want to run.")
        ("N_samples,Ns", po::value(&N_samples), "Number of samples taken during the simulation.")
        ("N_thermalization,Nth", po::value(&N_thermalization), "Number of thermalization steps.")
        ("N_between_samples,Nbs", po::value(&N_between_samples), "Number of steps between samples.")
        ("beta,bet", po::value(&beta), "Inverse temperature.")
        ("h_constant,hc", po::value(&h_constant), "Value of the constant h (electric field term).")
        ("h_hysteresis,hhys", po::value(&h_hys)->multitoken(), "Hysteresis values of the constant h (electric field term).")
        ("h_constant_therm,hct", po::value(&h_constant_therm), "Value of the constant h for thermalization (electric field term).")
        ("mu_constant,muc", po::value(&mu_constant), "Value of the constant mu_constant (star term).")
        ("J_constant,Jc", po::value(&J_constant), "Value of the constant J (plaquette term).")
        ("lmbda_constant,lmbdac", po::value(&lmbda_constant), "Value of the constant lmbda (gauge field term).")
        ("lmbda_hysteresis,lmbdahys", po::value(&lmbda_hys)->multitoken(), "Hysteresis values of the constant lmbda (gauge field term).")
        ("lmbda_constant_therm,lmbdact", po::value(&lmbda_constant_therm), "Value of the constant lmbda for thermalization (gauge field term).")
        ("N_resamples,Nr", po::value(&N_resamples), "Number of bootstrap resamples.")
        ("custom_therm,cth", po::value(&custom_therm), "Whether custom thermalization is used (to probe hysteresis).")
        ("observables,obs", po::value(&obs)->multitoken(), "Observables that are measured.") // can take multiple observables at once
        ("seed,s", po::value(&seed), "Seed for the pseudorandom number generator.")
        ("basis,bas", po::value(&basis), "Spin basis (\"x\" or \"z\").")
        ("output_directory,outdir", po::value(&path_output_directory), "Directory where the output is stored.")
        ("folder_name,fn", po::value(&folder_name), "Name of the output subfolder.")
        ("folder_names,fns", po::value(&folder_names)->multitoken(), "Directories where hysteresis output is stored.")
        ("snapshots,snap", po::value(&save_snapshots), "Whether snapshots should be saved.")
        ("full_time_series,fts", po::value(&full_time_series), "Whether full time series should be saved.")
        ("process_index,procid", po::value(&process_index), "Identifier of process (for debugging).")
        ("lattice_type,lat", po::value(&lattice_type), "Type of lattice used.")
        ("system_size,L", po::value(&system_size), "System size of lattice (one coordinate).")
        ("boundaries,bound", po::value(&boundaries), "Boundary condition of the lattice (periodic or open).")
        ("default_spin,dsp", po::value(&default_spin), "Default spin (electric field) for lattice initialization.");

        po::variables_map vm;
        po::command_line_parser parser(argc, argv);
        parser.options(desc).style(
            po::command_line_style::default_style
            | po::command_line_style::allow_long_disguise
        );
        po::store(parser.run(), vm);
        po::notify(vm); 

        BOOST_LOG_TRIVIAL(info) << "Subprocess [ID=" << process_index 
                                << "] SPAWNED with [β=" << beta 
                                << "], [μ=" << mu_constant 
                                << "], [J=" << J_constant 
                                << "], [h=" << h_constant 
                                << "], [λ=" << lmbda_constant << "].";

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        paratoric::LatSpec lat_spec = {
            basis,
            lattice_type,
            system_size,
            beta,
            boundaries,
            default_spin
        };
        paratoric::ParamSpec param_spec = {
            .mu = mu_constant,
            .h = h_constant,
            .J = J_constant,
            .lmbda = lmbda_constant,
            .h_therm = h_constant_therm,
            .lmbda_therm = lmbda_constant_therm,
            .h_hys = h_hys,
            .lmbda_hys = lmbda_hys
        };

        paratoric::SimSpec sim_spec = {
            .N_samples = N_samples,
            .N_thermalization = N_thermalization,
            .N_between_samples = N_between_samples,
            .N_resamples = N_resamples,
            .custom_therm = custom_therm,
            .seed = seed,
            .observables = obs
        };
        paratoric::OutSpec out_spec = {
            .path_out = path_output_directory,
            .folder_name = folder_name,
            .folder_names = folder_names,
            .save_snapshots = save_snapshots,
            .full_time_series = full_time_series
        };

        if (simulation == "etc_sample") {
            auto io = paratoric::IO();
            io.etc_sample(
                paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
        } else if (simulation == "etc_hysteresis") {
            auto io = paratoric::IO();
            io.etc_hysteresis(
                paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
        } else if (simulation == "etc_thermalization") {
            auto io = paratoric::IO();
            io.etc_thermalization(
                paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
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