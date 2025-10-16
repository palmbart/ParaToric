// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "paratoric/types/types.hpp"

#include <complex>
#include <filesystem>
#include <memory>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace paratoric {

class ExtendedToricCode {
public:
    ExtendedToricCode();
    ~ExtendedToricCode();

    ExtendedToricCode(ExtendedToricCode&&) noexcept;
    ExtendedToricCode& operator=(ExtendedToricCode&&) noexcept;

    ExtendedToricCode(const ExtendedToricCode&) = delete;
    ExtendedToricCode& operator=(const ExtendedToricCode&) = delete;

    /**
     * @brief This method will run a QMC thermalization of the extended toric code with the specified parameters and return observables and acceptance ratio diagnostics.
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
     * @param config.lat_spec.basis eigenbasis of the spins, either 'x' or 'z'
     * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
     * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
     * @param config.lat_spec.beta inverse temperature
     * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
     * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
     * @param config.out_spec.path_out output directory for snapshots 
     * @param config.out_spec.save_snapshots whether snapshots should be saved (every 10000th snapshot will be saved)
     * 
     * @return Result             the result object.
     * @return Result.series      Time series (per snapshot) of all requested observables,
     *                            measured during thermalization. Each observable’s entry
     *                            contains one value per recorded snapshot, in time order.
     * @return Result.acc_ratio   Time series of Monte Carlo acceptance ratios.
     * 
     * @throws std::invalid_argument  If an input is inconsistent (e.g., beta <= 0,
     *                                unknown basis, unsupported lattice/boundary).
     * @throws std::runtime_error     On RNG initialization failure or I/O errors
     *                                when saving snapshots.
     * 
     * @pre config.sim_spec.N_thermalization > 0
     * @pre config.lat_spec.beta > 0
     * @pre config.lat_spec.basis in {'x','z'}
     * 
     */
    static Result get_thermalization(const Config& config);

    /**
     * @brief This method will run a QMC simulation of the extended toric code with the specified parameters and return observables.
     * 
     * @param config the configuration object
     * @param config.sim_spec.N_samples the number of snapshots
     * @param config.sim_spec.N_thermalization the number of thermalization steps
     * @param config.sim_spec.N_between_samples the number of steps between snapshots
     * @param config.sim_spec.N_resamples the number of bootstrap resamples
     * @param config.sim_spec.custom_therm if custom thermalization is used (to probe hysteresis)
     * @param config.sim_spec.observables all observables that are calculated for every snapshot
     * @param config.sim_spec.seed the seed for the pseudorandom number generator
     * @param config.param_spec.mu the Hamiltonian parameter (star term)
     * @param config.param_spec.h the Hamiltonian parameter (electric field term)
     * @param config.param_spec.J the Hamiltonian parameter (plaquette term)
     * @param config.param_spec.lmbda the Hamiltonian parameter (gauge field term)
     * @param config.param_spec.h_therm the thermalization value of h, used if custom_therm enabled
     * @param config.param_spec.lmbda_therm the thermalization value of lmbda, used if custom_therm enabled
     * @param config.lat_spec.basis eigenbasis of the spins, either 'x' or 'z'
     * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
     * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
     * @param config.lat_spec.beta inverse temperature
     * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
     * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
     * @param config.out_spec.path_out output directory for snapshots 
     * @param config.out_spec.save_snapshots whether snapshots should be saved (every snapshot will be saved)
     * 
     * @return Result             the result object.
     * @return Result.series      Full counting statistics of all requested observables in the input order,
     *                            Each observable’s entry contains one value per recorded snapshot, 
     *                            in time order.
     * @return Result.mean        Bootstrap observable means
     * @return Result.mean_std    Bootstrap standard errors of the mean
     * @return Result.binder      Bootstrap binder ratios
     * @return Result.binder_std  Bootstrap standard errors of the binder ratios
     * @return Result.tau_int     Estimated integrated autocorrelation times
     * 
     * @throws std::invalid_argument  If an input is inconsistent (e.g., beta <= 0,
     *                                unknown basis, unsupported lattice/boundary).
     * @throws std::runtime_error     On RNG initialization failure or I/O errors
     *                                when saving snapshots.
     * 
     * @pre config.sim_spec.N_samples > 0
     * @pre config.sim_spec.N_thermalization > 0
     * @pre config.lat_spec.beta > 0
     * @pre config.lat_spec.basis in {'x','z'}
     * 
     */
    static Result get_sample(const Config& config);

    /**
     * @brief This method will run a QMC hysteresis simulation of the extended toric code with the specified parameters and return observables.
     * 
     * @param config the configuration object
     * @param config.sim_spec.N_samples the number of snapshots
     * @param config.sim_spec.N_thermalization the number of thermalization steps
     * @param config.sim_spec.N_between_samples the number of steps between snapshots
     * @param config.sim_spec.N_resamples the number of bootstrap resamples
     * @param config.sim_spec.observables all observables that are calculated for every snapshot
     * @param config.sim_spec.seed the seed for the pseudorandom number generator
     * @param config.param_spec.mu the Hamiltonian parameter (star term)
     * @param config.param_spec.h_hys the Hamiltonian parameters (electric field term)
     * @param config.param_spec.J the Hamiltonian parameter (plaquette term)
     * @param config.param_spec.lmbda_hys the Hamiltonian parameters (gauge field term)
     * @param config.lat_spec.basis eigenbasis of the spins, either 'x' or 'z'
     * @param config.lat_spec.lattice_type the lattice type, e.g. "triangular"
     * @param config.lat_spec.system_size the system size of the lattice (in one dimension)
     * @param config.lat_spec.beta inverse temperature
     * @param config.lat_spec.boundaries the boundary condition of the lattice (periodic, open)
     * @param config.lat_spec.default_spin the default spin on the links (1 or -1)
     * @param config.out_spec.paths_out all output directories for snapshots 
     * @param config.out_spec.save_snapshots whether snapshots should be saved (every snapshot will be saved)
     * 
     * @return Result                 the result object.
     * @return Result.series_hys      Full counting statistics of all requested observables.
     *                                Each element in the outer vector represents one parameter point (h_hys, lmbda_hys)
     *                                Each element in the middle vectors represent one observables in the order of input.
     *                                The inner vector are the full counting statistics in time order.
     * @return Result.mean_hys        Bootstrap observable means. 
     *                                Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
     *                                Each element in the inner vectors represents one observables in the order of input.
     * @return Result.mean_std_hys    Bootstrap standard errors of the mean.
     *                                Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
     *                                Each element in the inner vectors represents one observables in the order of input.
     * @return Result.binder_hys      Bootstrap binder ratios.
     *                                Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
     *                                Each element in the inner vectors represents one observables in the order of input.
     * @return Result.binder_std_hys  Bootstrap standard errors of the binder ratios.
     *                                Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
     *                                Each element in the inner vectors represents one observables in the order of input.
     * @return Result.tau_int_hys     Estimated integrated autocorrelation times.
     *                                Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
     *                                Each element in the inner vectors represents one observables in the order of input.
     * 
     * @throws std::invalid_argument  If an input is inconsistent (e.g., beta <= 0,
     *                                unknown basis, unsupported lattice/boundary).
     * @throws std::runtime_error     On RNG initialization failure or I/O errors
     *                                when saving snapshots.
     * 
     * @pre config.sim_spec.N_samples > 0
     * @pre config.sim_spec.N_thermalization > 0
     * @pre config.lat_spec.beta > 0
     * @pre config.lat_spec.basis in {'x','z'}
     * 
     */
    static Result get_hysteresis(const Config& config);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace paratoric
