// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "io/io.hpp"

#include "mcmc/extended_toric_code_qmc.hpp"
#include "paratoric/types/types.hpp"

#include <H5Cpp.h>  

#include <algorithm>   // std::min
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>  

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace paratoric {

void IO::etc_sample(
    const Config& config
    ) {

    std::string folder_name_new;
    if (!config.out_spec.folder_name.empty()){
        folder_name_new = config.out_spec.folder_name;
    } else {
        folder_name_new = config.lat_spec.lattice_type + "_" + std::to_string(config.lat_spec.system_size) + "_" + config.lat_spec.boundaries + "_"
            + std::to_string(config.lat_spec.beta);
    }
    std::filesystem::path path_folder_name(folder_name_new);
    std::filesystem::path path_out = config.out_spec.path_out / path_folder_name;
    std::filesystem::create_directories(path_out);
    
    Result result_spec;
    std::vector<std::string> obs_types;
    if (config.lat_spec.basis == 'x') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'x'>>();
        result_spec = mc->get_sample(
            Config{config.sim_spec, config.param_spec, config.lat_spec, 
                OutSpec{.path_out=path_out, .save_snapshots=config.out_spec.save_snapshots}}
        ); 
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    } else if (config.lat_spec.basis == 'z') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'z'>>();
        result_spec = mc->get_sample(
            Config{config.sim_spec, config.param_spec, config.lat_spec, 
                OutSpec{.path_out=path_out, .save_snapshots=config.out_spec.save_snapshots}}
        ); 
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    }
    const auto& result = result_spec.series;
    const auto& obs_means = result_spec.mean;
    const auto& obs_std = result_spec.mean_std;
    const auto& binders_means = result_spec.binder;
    const auto& binders_std = result_spec.binder_std;
    const auto& autocorrelation_time = result_spec.tau_int;

    std::filesystem::path hdf5_name("obs.h5");
    std::filesystem::path hdf5_path_out = path_out / hdf5_name;

    H5::H5File file{ hdf5_path_out, H5F_ACC_TRUNC };

    H5::Group sim_grp = getOrCreateGroup(file, "simulation");
    H5::Group results_grp = getOrCreateGroup(sim_grp, "results");

    for (size_t i = 0; i < config.sim_spec.observables.size(); ++i) {
        auto obs_name = config.sim_spec.observables[i];
        auto obs_type = obs_types[i];
        auto obs_result_vector = result[i];

        H5::Group obs_grp = getOrCreateGroup(results_grp, obs_name);

        if (obs_type == "real" || obs_type == "fredenhagen_marcu" || obs_type == "susceptibility") {
            if (config.out_spec.full_time_series) {
                hsize_t dims[1] = { obs_result_vector.size() };
                H5::DataSpace dataspace{ 1, dims };

                H5::CompType complex_data_type(sizeof(obs_result_vector[0]));
                complex_data_type.insertMember( "r", 0, H5::PredType::NATIVE_DOUBLE);
                complex_data_type.insertMember( "i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

                H5::DataSet dataset = obs_grp.createDataSet("series", complex_data_type, dataspace);

                dataset.write(obs_result_vector.data(), complex_data_type);
            }
        } else {
            throw std::runtime_error(std::format("Observable type \"{}\" is not supported.", obs_type));
        }

        H5::DataSpace scalar_space(H5S_SCALAR);
        auto ds_means = obs_grp.createDataSet("mean", H5::PredType::NATIVE_DOUBLE, scalar_space);
        ds_means.write(&(obs_means[i]), H5::PredType::NATIVE_DOUBLE);
        auto ds_mean_std = obs_grp.createDataSet("mean_error", H5::PredType::NATIVE_DOUBLE, scalar_space);
        ds_mean_std.write(&(obs_std[i]), H5::PredType::NATIVE_DOUBLE);
        auto ds_binder = obs_grp.createDataSet("binder", H5::PredType::NATIVE_DOUBLE, scalar_space);
        ds_binder.write(&(binders_means[i]), H5::PredType::NATIVE_DOUBLE);
        auto ds_binder_std = obs_grp.createDataSet("binder_error", H5::PredType::NATIVE_DOUBLE, scalar_space);
        ds_binder_std.write(&(binders_std[i]), H5::PredType::NATIVE_DOUBLE);
        auto ds_autocorrelation = obs_grp.createDataSet("autocorrelation_time", H5::PredType::NATIVE_DOUBLE, scalar_space);
        ds_autocorrelation.write(&(autocorrelation_time[i]), H5::PredType::NATIVE_DOUBLE);
    }
}

void IO::etc_hysteresis(
    const Config& config
    ) {

    std::vector<std::filesystem::path> paths_out;

    for (auto& folder_name : config.out_spec.folder_names) {
        std::filesystem::path path_folder_name(folder_name);
        std::filesystem::path path_out = config.out_spec.path_out / path_folder_name;
        std::filesystem::create_directories(path_out);
        paths_out.emplace_back( std::move(path_out) );
    }

    Result result_spec;
    std::vector<std::string> obs_types;
    if (config.lat_spec.basis == 'x') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'x'>>();
        result_spec = mc->get_hysteresis(
            Config{config.sim_spec, config.param_spec, config.lat_spec, OutSpec{.paths_out=paths_out, .save_snapshots=config.out_spec.save_snapshots}}
        );
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    } else if (config.lat_spec.basis == 'z') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'z'>>();
        result_spec = mc->get_hysteresis(
            Config{config.sim_spec, config.param_spec, config.lat_spec, OutSpec{.paths_out=paths_out, .save_snapshots=config.out_spec.save_snapshots}}
        );
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    }

    const auto& results = result_spec.series_hys;
    const auto& hys_means = result_spec.mean_hys;
    const auto& hys_means_std = result_spec.mean_std_hys;
    const auto& hys_binders = result_spec.binder_hys;
    const auto& hys_binders_std = result_spec.binder_std_hys;
    const auto& hys_autocorrelation_time = result_spec.tau_int_hys;

    for (size_t n = 0; n < std::min( config.param_spec.h_hys.size(), std::min(config.param_spec.lmbda_hys.size(), paths_out.size()) ); n++) {
        auto path_out = paths_out[ n ];
        auto result = results[ n ];
        auto obs_means = hys_means[ n ];
        auto obs_std = hys_means_std[ n ];
        auto binders_means = hys_binders[ n ];
        auto binders_std = hys_binders_std[ n ];
        auto autocorrelation_time = hys_autocorrelation_time[ n ];

        std::filesystem::path hdf5_name("obs.h5");
        std::filesystem::path hdf5_path_out = path_out / hdf5_name;

        H5::H5File file{ hdf5_path_out, H5F_ACC_TRUNC };

        H5::Group sim_grp = getOrCreateGroup(file, "simulation");
        H5::Group results_grp = getOrCreateGroup(sim_grp, "results");

        for (size_t i = 0; i < config.sim_spec.observables.size(); ++i) {
            auto obs_name = config.sim_spec.observables[i];
            auto obs_type = obs_types[i];
            auto obs_result_vector = result[i];

            H5::Group obs_grp = getOrCreateGroup(results_grp, obs_name);

            if (obs_type == "real" || obs_type == "fredenhagen_marcu" || obs_type == "susceptibility") {
                if (config.out_spec.full_time_series) {
                    hsize_t dims[1] = { obs_result_vector.size() };
                    H5::DataSpace dataspace{ 1, dims };

                    H5::CompType complex_data_type(sizeof(obs_result_vector[0]));
                    complex_data_type.insertMember( "r", 0, H5::PredType::NATIVE_DOUBLE);
                    complex_data_type.insertMember( "i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

                    H5::DataSet dataset = obs_grp.createDataSet("series", complex_data_type, dataspace);

                    dataset.write(obs_result_vector.data(), complex_data_type);
                }
            } else {
                throw std::runtime_error(std::format("Observable type \"{}\" is not supported.", obs_type));
            }

            H5::DataSpace scalar_space(H5S_SCALAR);
            auto ds_means = obs_grp.createDataSet("mean", H5::PredType::NATIVE_DOUBLE, scalar_space);
            ds_means.write(&(obs_means[i]), H5::PredType::NATIVE_DOUBLE);
            auto ds_mean_std = obs_grp.createDataSet("mean_error", H5::PredType::NATIVE_DOUBLE, scalar_space);
            ds_mean_std.write(&(obs_std[i]), H5::PredType::NATIVE_DOUBLE);
            auto ds_binder = obs_grp.createDataSet("binder", H5::PredType::NATIVE_DOUBLE, scalar_space);
            ds_binder.write(&(binders_means[i]), H5::PredType::NATIVE_DOUBLE);
            auto ds_binder_std = obs_grp.createDataSet("binder_error", H5::PredType::NATIVE_DOUBLE, scalar_space);
            ds_binder_std.write(&(binders_std[i]), H5::PredType::NATIVE_DOUBLE);
            auto ds_autocorrelation_time = obs_grp.createDataSet("autocorrelation_time", H5::PredType::NATIVE_DOUBLE, scalar_space);
            ds_autocorrelation_time.write(&(autocorrelation_time[i]), H5::PredType::NATIVE_DOUBLE);
        }
    }
}

void IO::etc_thermalization(
    const Config& config
    ) { 

    std::string folder_name_new;
    if (!config.out_spec.folder_name.empty()){
        folder_name_new = config.out_spec.folder_name;
    } else {
        folder_name_new = config.lat_spec.lattice_type + "_" + std::to_string(config.lat_spec.system_size) + "_" + config.lat_spec.boundaries + "_"
            + std::to_string(config.lat_spec.beta);
    }
    std::filesystem::path path_folder_name(folder_name_new);
    std::filesystem::path path_out = config.out_spec.path_out / path_folder_name;
    std::filesystem::create_directories(path_out);

    Result result_spec;
    std::vector<std::string> obs_types;
    if (config.lat_spec.basis == 'x') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'x'>>();
        result_spec = result_spec = mc->get_thermalization(
            Config{config.sim_spec, config.param_spec, config.lat_spec, OutSpec{.path_out=path_out, .save_snapshots=config.out_spec.save_snapshots}}
        );
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    } else if (config.lat_spec.basis == 'z') {
        auto mc = std::make_unique<ExtendedToricCodeQMC<'z'>>();
        result_spec = result_spec = mc->get_thermalization(
            Config{config.sim_spec, config.param_spec, config.lat_spec, OutSpec{.path_out=path_out, .save_snapshots=config.out_spec.save_snapshots}}
        );
        obs_types = mc->get_obs_type_vec(config.sim_spec.observables);
    }

    const auto& obs_result = result_spec.series;
    const auto& acc_ratio_result = result_spec.acc_ratio;

    std::filesystem::path hdf5_name("obs.h5");
    std::filesystem::path hdf5_path_out = path_out / hdf5_name;

    H5::H5File file{ hdf5_path_out, H5F_ACC_TRUNC };

    H5::Group sim_grp = getOrCreateGroup(file, "simulation");
    H5::Group results_grp = getOrCreateGroup(sim_grp, "results");

    for (size_t i = 0; i < config.sim_spec.observables.size(); ++i) {
        auto obs_name = config.sim_spec.observables[i];
        auto obs_type = obs_types[i];
        auto obs_result_vector = obs_result[i];

        H5::Group obs_grp = getOrCreateGroup(results_grp, obs_name);

        if (obs_type == "real" || obs_type == "fredenhagen_marcu" || obs_type == "susceptibility") {
            hsize_t dims[1] = { obs_result_vector.size() };
            H5::DataSpace dataspace{ 1, dims };

            H5::CompType complex_data_type(sizeof(obs_result_vector[0]));
            complex_data_type.insertMember( "r", 0, H5::PredType::NATIVE_DOUBLE);
            complex_data_type.insertMember( "i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

            H5::DataSet dataset = obs_grp.createDataSet("series", complex_data_type, dataspace);

            dataset.write(obs_result_vector.data(), complex_data_type);
        } else {
            throw std::runtime_error(std::format("Observable type \"{}\" is not supported.", obs_type));
        }
    }


    hsize_t dims[1] = { acc_ratio_result.size() };
    H5::DataSpace dataspace{ 1, dims };

    auto dtype = H5::PredType::NATIVE_DOUBLE;

    H5::DataSet dataset = results_grp.createDataSet("acc_ratio", dtype, dataspace);

    dataset.write(acc_ratio_result.data(), dtype);
}

} // namespace paratoric