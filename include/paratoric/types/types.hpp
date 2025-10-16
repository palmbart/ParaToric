// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <complex>
#include <filesystem>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace paratoric {

/**
 * @brief Lattice specification
 * 
 * @param basis eigenbasis of the spins, either 'x' or 'z'
 * @param lattice_type the lattice type, e.g. "triangular"
 * @param system_size the system size of the lattice (in one dimension)
 * @param beta inverse temperature
 * @param boundaries the boundary condition of the lattice (periodic, open)
 * @param default_spin the default spin on the links (1 or -1)
 * 
 */
struct LatSpec {
    char        basis        = 'x';
    std::string lattice_type = "square";
    int         system_size  = 16;
    double      beta         = 16.;
    std::string boundaries   = "periodic";
    int         default_spin = +1;
};

/**
 * @brief Hamiltonian parameter specification
 * 
 * @param mu the Hamiltonian parameter (star term)
 * @param h the Hamiltonian parameter (electric field term)
 * @param J the Hamiltonian parameter (plaquette term)
 * @param lmbda the Hamiltonian parameter (gauge field term)
 * @param h_therm the thermalization value of h, used if custom_therm enabled
 * @param lmbda_therm the thermalization value of lmbda, used if custom_therm enabled
 * @param h_hys the Hamiltonian hysteresis parameters (electric field term)
 * @param lmbda_hys the Hamiltonian hysteresis parameters (gauge field term)
 * 
 */
struct ParamSpec {
    double mu     = 1.0;
    double h      = 0.0;
    double J      = 1.0;
    double lmbda  = 0.0;
    // optional “thermalization schedule” values when they differ
    double h_therm     = std::numeric_limits<double>::quiet_NaN();
    double lmbda_therm = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> h_hys{};
    std::vector<double> lmbda_hys{};
};

/**
 * @brief Simulation specification
 * 
 * @param N_samples the number of snapshots
 * @param N_thermalization the number of thermalization steps
 * @param N_between_samples the number of steps between snapshots
 * @param N_resamples the number of bootstrap resamples
 * @param custom_therm if custom thermalization is used (to probe hysteresis)
 * @param seed the seed for the pseudorandom number generator
 * @param observables all observables that are calculated for every snapshot
 * 
 */
struct SimSpec {
    int    N_samples                = 1000;     // 0 for thermalization-only
    int    N_thermalization         = 10000;
    int    N_between_samples        = 1000;
    int    N_resamples              = 1000;
    bool   custom_therm             = false; // used only by sampling
    int    seed                     = 0;     // 0 means random seed
    std::vector<std::string> observables{};
};

/**
 * @brief Output specification
 * 
 * @param path_out the directory where output files will be written to 
 * @param paths_out the directories (hysteresis simulation) where output files will be written to
 * @param folder_name the folder name of the output folder, additional possibility to make the output clearer
 * @param folder_names the folder names of the output folders (hysteresis simulation), additional possibility to make the output clearer
 * @param save_snapshots whether snapshots should be saved
 * @param full_time_series whether full time series should be saved
 * 
 */
struct OutSpec {
    std::filesystem::path path_out; // single path
    std::vector<std::filesystem::path> paths_out{}; // hysteresis
    std::string folder_name; // single folder
    std::vector<std::string> folder_names{}; // hysteresis
    bool save_snapshots = false;
    bool full_time_series = 0;
};

/**
 * @brief Configuration for starting a simulation
 * 
 * @param sim_spec Simulation parameters (N_samples, seed, ...) 
 * @param param_spec Hamiltonian parameters
 * @param lat_spec Lattice parameters (system_size, beta, ...)
 * @param out_spec Output parameters (path_out, folder_name, ...)
 * 
 */
struct Config{
    SimSpec sim_spec{};
    ParamSpec param_spec{};
    LatSpec lat_spec{};
    OutSpec out_spec{};
};

/**
 * @brief Result of a simulation
 * 
 * @param series          Full counting statistics of all requested observables in the input order,
 *                        Each observable’s entry contains one value per recorded snapshot, 
 *                        in time order.
 * @param acc_ratio       Time series of Monte Carlo acceptance ratios.
 * @param mean            Bootstrap observable means
 * @param mean_std        Bootstrap standard errors of the mean
 * @param binder          Bootstrap binder ratios
 * @param binder_std      Bootstrap standard errors of the binder ratios
 * @param tau_int         Estimated integrated autocorrelation times
 * @param series_hys      Full counting statistics of all requested observables.
 *                        Each element in the outer vector represents one parameter point (h_hys, lmbda_hys)
 *                        Each element in the middle vectors represent one observables in the order of input.
 *                        The inner vector are the full counting statistics in time order.
 * @param mean_hys        Bootstrap observable means. 
 *                        Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                        Each element in the inner vectors represents one observables in the order of input.
 * @param mean_std_hys    Bootstrap standard errors of the mean.
 *                        Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                        Each element in the inner vectors represents one observables in the order of input.
 * @param binder_hys      Bootstrap binder ratios.
 *                        Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                        Each element in the inner vectors represents one observables in the order of input.
 * @param binder_std_hys  Bootstrap standard errors of the binder ratios.
 *                        Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                        Each element in the inner vectors represents one observables in the order of input.
 * @param tau_int_hys     Estimated integrated autocorrelation times.
 *                        Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                        Each element in the inner vectors represents one observables in the order of input.
 * 
 */
struct Result {
    std::vector<std::vector<std::variant<std::complex<double>, double>>> series{};
    std::vector<double> acc_ratio{};
    std::vector<double> mean{}, mean_std{}, binder{}, binder_std{}, tau_int{};
    std::vector<std::vector<std::vector<std::variant<std::complex<double>, double>>>> series_hys{};
    std::vector<std::vector<double>> mean_hys{}, mean_std_hys{}, binder_hys{}, binder_std_hys{}, tau_int_hys{};
};

} // namespace paratoric