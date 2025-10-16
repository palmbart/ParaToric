// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "paratoric/types/types.hpp"

#include <H5Cpp.h>   

#include <filesystem>
#include <vector> 

namespace paratoric {

class IO {
    public:
        IO() = default;
        ~IO() = default;

        IO(IO const&) = default;
        IO& operator=(IO const&) = default;

        /**
         * @brief This method will run a QMC simulation of the extended toric code with the specified parameters and write the results to a HDF5 file. The Thermalization results are discarded.
         * 
         * @param config the configuration object
         * @param config.sim_spec.N_samples the number of snapshots
         * @param config.sim_spec.N_thermalization the number of thermalization steps
         * @param config.sim_spec.N_between_samples the number of steps between snapshots
         * @param config.sim_spec.N_resamples the number of bootstrap resamples
         * @param config.sim_spec.custom_therm whether custom thermalization is used (to probe hysteresis)
         * @param config.sim_spec.observables all observables that are calculated for every snapshot
         * @param config.sim_spec.seed the seed for the pseudorandom number generator
         * @param config.param_spec.mu the Hamiltonian parameter (star term)
         * @param config.param_spec.h the Hamiltonian parameter (electric field term)
         * @param config.param_spec.J the Hamiltonian parameter (plaquette term)
         * @param config.param_spec.lmbda the Hamiltonian parameter (gauge field term)
         * @param config.lat_spec.basis spin basis (x or z)
         * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
         * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
         * @param config.lat_spec.beta inverse temperature
         * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
         * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
         * @param config.out_spec.path_out the directory where the method will write the output files to 
         * @param config.out_spec.folder_name the folder name of the output folder, additional possibility to make the output clearer
         * @param config.out_spec.save_snapshots whether snapshots should be saved
         * @param config.out_spec.full_time_series whether full time series should be saved
         * 
         */
        static void etc_sample(
            const Config& config
        ); 

        /**
         * @brief This method will run a QMC hysteresis simulation of the extended toric code with the specified parameters and write the results to a HDF5 file. The Thermalization results are discarded.
         * 
         * @param config the configuration object
         * @param config.sim_spec.N_samples the number of snapshots
         * @param config.sim_spec.N_thermalization the number of thermalization steps
         * @param config.sim_spec.N_between_snapshots the number of samples between snapshots
         * @param config.sim_spec.N_resamples the number of bootstrap resamples
         * @param config.sim_spec.observables all observables that are calculated for every snapshot
         * @param config.sim_spec.seed the seed for the pseudorandom number generator
         * @param config.param_spec.mu the Hamiltonian parameter (star term)
         * @param config.param_spec.h_hys the Hamiltonian parameters (electric field term)
         * @param config.param_spec.J the Hamiltonian parameter (plaquette term)
         * @param config.param_spec.lmbda_hys the Hamiltonian parameters (gauge field term)
         * @param config.lat_spec.basis spin basis (x or z)
         * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
         * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
         * @param config.lat_spec.beta inverse temperature
         * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
         * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
         * @param config.out_spec.paths_out all the directories where the method will write the output files to 
         * @param config.out_spec.folder_names the folder names of the output folder, additional possibility to make the output clearer
         * @param config.out_spec.save_snapshots whether snapshots should be saved
         * @param config.out_spec.full_time_series whether full time series should be saved
         * 
         */
        static void etc_hysteresis(
            const Config& config
        ); 

        /**
         * @brief This method will run a QMC simulation of the extended toric code with the specified parameters and write the results to a HDF5 file. Here we specifically look at the thermalization phase.
         * 
         * @param config the configuration object
         * @param config.sim_spec.N_thermalization the number of thermalization steps
         * @param config.sim_spec.N_resamples the number of bootstrap resamples
         * @param config.sim_spec.observables all observables that are calculated for every snapshot
         * @param config.sim_spec.seed the seed for the pseudorandom number generator
         * @param config.param_spec.mu the Hamiltonian parameter (star term)
         * @param config.param_spec.h the Hamiltonian parameter (electric field term)
         * @param config.param_spec.J the Hamiltonian parameter (plaquette term)
         * @param config.param_spec.lmbda the Hamiltonian parameter (gauge field term)
         * @param config.lat_spec.basis spin basis (x or z)
         * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
         * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
         * @param config.lat_spec.beta inverse temperature
         * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
         * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
         * @param config.out_spec.path_out the directory where the method will write the output files to 
         * @param config.out_spec.folder_name the folder name of the output folder, additional possibility to make the output clearer
         * @param config.out_spec.save_snapshots whether snapshots should be saved (every 10000th snapshot will be saved)
         * 
         */
        static void etc_thermalization(
            const Config& config
        );

    private:
        template<typename Parent>
        static H5::Group getOrCreateGroup(Parent& parent, const std::string& name) {
            htri_t exists = H5Lexists(parent.getLocId(), name.c_str(), H5P_DEFAULT);
            if (exists > 0) {
                return parent.openGroup(name);
            }
            // either not existing or an error, just create
            return parent.createGroup(name);
        }
};

} // namespace paratoric