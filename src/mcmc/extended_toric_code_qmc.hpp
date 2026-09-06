// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2026  Simon Mathias Linsel, Lode Pollet

#pragma once

#include "lattice/lattice.hpp"
#include "mcmc/input_validation.hpp"
#include "paratoric/types/types.hpp"
#include "rng/rng.hpp"
#include "statistics/autocorrelation.hpp"
#include "statistics/bootstrap.hpp"

#include <boost/container/small_vector.hpp>
#include <boost/log/core.hpp> 
#include <boost/log/expressions.hpp> 
#include <boost/log/trivial.hpp> 

#include <algorithm> 
#include <chrono>
#include <cmath>
#include <concepts>
#include <complex>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <span>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#define UNUSED(expr) do { (void)(expr); } while (0)

namespace paratoric {

/// Supported compile-time spin bases.
template<char B>
concept ValidBasis = (B == 'x' || B == 'z');

/**
 * @brief Continuous-time Metropolis sampler for the extended toric code.
 * @tparam Basis Spin eigenbasis, 'x' or 'z'; must match Config::lat_spec.basis.
 *
 * In the x-basis, edge and star terms are diagonal; single-spin and plaquette
 * flips carry lmbda and J. In the z-basis, edge and plaquette terms are diagonal;
 * single-spin and star flips carry h and mu. An update's "tuple" is a plaquette
 * in the x-basis and a star in the z-basis.
 *
 * Lattice stores bare integrals of spins and spin products. This class applies
 * the Hamiltonian minus signs and couplings to form the integrated potential
 * energy and Metropolis ratios. Accepted proposals update event histories,
 * the affected active caches, and the running energy together. A zero coupling
 * may leave its bare cache stale; rebuild caches before changing couplings.
 *
 * Each workflow constructs its own lattice and shares this backend's RNG with
 * it and the bootstrap routines. Copying the backend shares the RNG as well.
 * A nonzero configured seed reseeds it; zero preserves its current state.
 */
template<char Basis>
requires ValidBasis<Basis>
class ExtendedToricCodeQMC {
    public:
        using RNG = paratoric::rng::RNG;
        using SmallIndexVector = Lattice::SmallIndexVector;
        using SmallEnergyVector = Lattice::SmallEnergyVector;
        using SmallBoolVector = boost::container::small_vector<bool, 8>;

        /** @brief Use the supplied RNG, or create one seeded by std::random_device. */
        ExtendedToricCodeQMC(std::shared_ptr<RNG> rng = nullptr) 
        : rng(rng ? std::move(rng) : std::make_shared<RNG>()) {};

        ~ExtendedToricCodeQMC() = default;

        ExtendedToricCodeQMC(ExtendedToricCodeQMC const&) = default;
        ExtendedToricCodeQMC& operator=(ExtendedToricCodeQMC const&) = default;

        /**
         * @brief This method will run a QMC thermalization of the extended toric code with the specified parameters and return observables and acceptance ratio diagnostics.
         * 
         * @tparam Basis eigenbasis of the spins, either 'x' or 'z'
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
         * @pre config.sim_spec.N_thermalization >= 0
         * @pre config.sim_spec.N_resamples > 0
         * @pre config.lat_spec.beta > 0
         * @pre config.lat_spec.basis in {'x','z'}
         * 
         */
        Result get_thermalization(
            const Config& config
        );

        /**
         * @brief This method will run a QMC simulation of the extended toric code with the specified parameters and return observables.
         * 
         * @tparam Basis eigenbasis of the spins, either 'x' or 'z'
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
         * @pre config.sim_spec.N_thermalization >= 0
         * @pre config.sim_spec.N_between_samples >= 0
         * @pre config.sim_spec.N_resamples > 0
         * @pre config.lat_spec.beta > 0
         * @pre config.lat_spec.basis in {'x','z'}
         * 
         */
        Result get_sample(
            const Config& config
        );

        /**
         * @brief This method will run a QMC hysteresis simulation of the extended toric code with the specified parameters and return observables.
         * 
         * @tparam Basis eigenbasis of the spins, either 'x' or 'z'
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
         * @pre config.sim_spec.N_thermalization >= 0
         * @pre config.sim_spec.N_between_samples >= 0
         * @pre config.sim_spec.N_resamples > 0
         * @pre config.param_spec.h_hys and config.param_spec.lmbda_hys are non-empty
         *      and have equal lengths
         * @pre config.lat_spec.beta > 0
         * @pre config.lat_spec.basis in {'x','z'}
         * 
         */
        Result get_hysteresis(
            const Config& config
        );
        
        /**
         * @brief Resolve observable statistics categories in input order.
         * @param observables Registered observable names.
         * @return Categories: "real", "fredenhagen_marcu", or "susceptibility".
         * @throws std::invalid_argument If any name is not registered.
         */
        std::vector<std::string> get_obs_type_vec(const std::vector<std::string>& observables);
    
    private:
        // Measurements receive (lattice, h, lmbda, mu, J). Complex results pack
        // paired real estimators for the statistics routines below.
        std::function<double(Lattice&, double, double, double, double)> 
        percolation_probability_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.percolation_probability(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        plaquette_percolation_probability_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.plaquette_percolation_probability(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        cube_percolation_probability 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.cube_percolation_probability(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        percolation_strength_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.percolation_strength(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        plaquette_percolation_strength_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.plaquette_percolation_strength(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        string_number_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.get_string_count(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        largest_cluster_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.largest_cluster(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        largest_plaquette_cluster_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.largest_plaquette_cluster(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        anyon_count_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.get_anyon_count(); 
        };

        std::function<std::complex<double>(Lattice&, double, double, double, double)> 
        fredenhagen_marcu_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            return lat.fredenhagen_marcu(); 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        staggered_imaginary_times_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') return lat.get_staggered_imaginary_times_plaquette(); 
            else return lat.get_staggered_imaginary_times_star();
        };

        std::function<double(Lattice&, double, double, double, double)> 
        energy_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                return - lat.get_diag_single_energy() * h - lat.get_diag_tuple_energy_x() * mu 
                - lat.get_non_diag_single_energy_x() - lat.get_non_diag_tuple_energy_x();
            }
            else {
                return - lat.get_diag_single_energy() * lmbda - lat.get_diag_tuple_energy_z() * J 
                - lat.get_non_diag_single_energy_z() - lat.get_non_diag_tuple_energy_z();
            } 
        };

        std::function<double(Lattice&, double, double, double, double)> 
        energy_h_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') return - lat.get_diag_single_energy() * h; 
            else return - lat.get_non_diag_single_energy_z();
        };

        std::function<double(Lattice&, double, double, double, double)> 
        energy_mu_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') return - lat.get_diag_tuple_energy_x() * mu; 
            else return - lat.get_non_diag_tuple_energy_z();
        };

        std::function<double(Lattice&, double, double, double, double)> 
        energy_lmbda_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') return - lat.get_non_diag_single_energy_x(); 
            else return - lat.get_diag_single_energy() * lmbda;
        };

        std::function<double(Lattice&, double, double, double, double)> 
        energy_J_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') return - lat.get_non_diag_tuple_energy_x(); 
            else return - lat.get_diag_tuple_energy_z() * J;
        }; 

        std::function<double(Lattice&, double, double, double, double)> 
        sigma_x_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                return lat.get_diag_single_energy()/static_cast<double>(lat.get_edge_count()); 
            } else {
                if (h == 0.) return 0.;
                return lat.get_non_diag_single_energy_z()/static_cast<double>(lat.get_edge_count() * h);
            }
        };

        std::function<std::complex<double>(Lattice&, double, double, double, double)> 
        sigma_x_static_susceptibility_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            lat.rotate_imag_time();
            if constexpr (Basis == 'x') return lat.get_diag_M_M(); 
            else return lat.get_non_diag_M_M();
        };

        std::function<std::complex<double>(Lattice&, double, double, double, double)> 
        sigma_x_dynamical_susceptibility_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            lat.rotate_imag_time();
            if constexpr (Basis == 'x') {
                return lat.get_diag_dynamical_M_M();
            } else {
                if (h == 0.) return std::complex<double>{};
                return lat.get_kL_kR_single() / static_cast<double>(std::sqrt(2) * h);
            }
        };

        std::function<double(Lattice&, double, double, double, double)> 
        sigma_z_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                if (lmbda == 0.) return 0.;
                return lat.get_non_diag_single_energy_x()/static_cast<double>(lat.get_edge_count() * lmbda); 
            } else {
                return (lat.get_diag_single_energy()/static_cast<double>(lat.get_edge_count()));
            }
        };

        std::function<std::complex<double>(Lattice&, double, double, double, double)> 
        sigma_z_static_susceptibility_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            lat.rotate_imag_time();
            if constexpr (Basis == 'x') return lat.get_non_diag_M_M(); 
            else return lat.get_diag_M_M();
        };

        std::function<std::complex<double>(Lattice&, double, double, double, double)> 
        sigma_z_dynamical_susceptibility_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            lat.rotate_imag_time();
            if constexpr (Basis == 'x') {
                if (lmbda == 0.) return std::complex<double>{};
                return lat.get_kL_kR_single() / static_cast<double>(std::sqrt(2) * lmbda);
            } else {
                return lat.get_diag_dynamical_M_M();
            }
        };

        std::function<double(Lattice&, double, double, double, double)> 
        star_x_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                return lat.get_diag_tuple_energy_x()/static_cast<double>(lat.get_vertex_count()); 
            } else {
                if (mu == 0.) return 0.;
                return lat.get_non_diag_tuple_energy_z()/static_cast<double>(lat.get_vertex_count() * mu);
            }
        };

        std::function<double(Lattice&, double, double, double, double)> 
        plaquette_z_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                if (J == 0.) return 0.;
                return lat.get_non_diag_tuple_energy_x()/static_cast<double>(lat.get_plaquette_count() * J); 
            } else {
                return lat.get_diag_tuple_energy_z()/static_cast<double>(lat.get_plaquette_count());
            }
        };

        std::function<double(Lattice&, double, double, double, double)> 
        delta_obs = [](Lattice& lat, double h, double lmbda, double mu, double J) { 
            if constexpr (Basis == 'x') {
                const double star_x = lat.get_diag_tuple_energy_x()
                    / static_cast<double>(lat.get_vertex_count());
                if (J == 0.) return -star_x;
                const double plaquette_z = lat.get_non_diag_tuple_energy_x()
                    / static_cast<double>(lat.get_plaquette_count() * J);
                return plaquette_z - star_x;
            } else {
                const double plaquette_z = lat.get_diag_tuple_energy_z()
                    / static_cast<double>(lat.get_plaquette_count());
                if (mu == 0.) return plaquette_z;
                const double star_x = lat.get_non_diag_tuple_energy_z()
                    / static_cast<double>(lat.get_vertex_count() * mu);
                return plaquette_z - star_x;
            }
        };

        std::function<double(Lattice&, double, double, double, double)> 
        anyon_density_obs 
        = [](Lattice& lat, double h, double lmbda, double mu, double J) {
            if constexpr (Basis == 'x') 
                return lat.get_anyon_count()/static_cast<double>(lat.get_vertex_count()); 
            else 
                return lat.get_anyon_count()/static_cast<double>(lat.get_plaquette_count()); 
        };

        /** @brief Registry entry connecting an observable to its estimator and statistics category. */
        struct Obs {
            public:
                std::string obs_name;
                std::string obs_type;
                std::function<
                    std::variant< std::complex<double>, double>(Lattice&, double, double, double, double)
                > obs_func;
        };

        // Keep each statistics category consistent with its estimator's packed values.
        std::vector<Obs> obs_vec = {
            {"anyon_count", "real", anyon_count_obs},
            {"anyon_density", "real", anyon_density_obs},
            {"cube_percolation_probability", "real", cube_percolation_probability},
            {"delta", "real", delta_obs},
            {"energy", "real", energy_obs},
            {"energy_h", "real", energy_h_obs},
            {"energy_lmbda", "real", energy_lmbda_obs},
            {"energy_J", "real", energy_J_obs},
            {"energy_mu", "real", energy_mu_obs},
            {"fredenhagen_marcu", "fredenhagen_marcu", fredenhagen_marcu_obs},
            {"largest_cluster", "real", largest_cluster_obs},
            {"largest_plaquette_cluster", "real", largest_plaquette_cluster_obs},
            {"percolation_probability", "real", percolation_probability_obs},
            {"percolation_strength", "real", percolation_strength_obs},
            {"plaquette_percolation_probability", "real", plaquette_percolation_probability_obs},
            {"plaquette_percolation_strength", "real", plaquette_percolation_strength_obs},
            {"plaquette_z", "real", plaquette_z_obs},
            {"sigma_x", "real", sigma_x_obs},
            {"sigma_x_static_susceptibility", "susceptibility", sigma_x_static_susceptibility_obs},
            {"sigma_x_dynamical_susceptibility", "susceptibility", sigma_x_dynamical_susceptibility_obs},
            {"sigma_z", "real", sigma_z_obs},
            {"sigma_z_static_susceptibility", "susceptibility", sigma_z_static_susceptibility_obs},
            {"sigma_z_dynamical_susceptibility", "susceptibility", sigma_z_dynamical_susceptibility_obs},
            {"staggered_imaginary_times", "real", staggered_imaginary_times_obs},
            {"star_x", "real", star_x_obs},
            {"string_number", "real", string_number_obs}
            };
        
        /**
         * @brief Resolve measurement functions in the requested observable order.
         * @param observables Registered observable names.
         * @return Functions taking (lattice, h, lmbda, mu, J).
         * @throws std::invalid_argument If any name is not registered.
         */
        std::vector<std::function<std::variant< std::complex<double>, double>(Lattice&, double, double, double, double)>>
        get_obs_func_vec(
            const std::vector<std::string>& observables
        );

        // Shared with the lattice and bootstrap so a configured seed covers the whole run.
        std::shared_ptr<RNG> rng;
        std::uniform_real_distribution<double> uniform_dist{0., 1.};
        static constexpr double PRECISION = std::numeric_limits<double>::epsilon();
        static constexpr double AUTOCORRELATION_WARNING_SAMPLE_FRACTION = 0.1;

        /** @brief Draw an index in [0, bound); bound must be positive. */
        int random_index(int bound) {
            return static_cast<int>(
                paratoric::rng::uniform_index(*rng, static_cast<std::uint64_t>(bound))
            );
        }

        double uniform_real(double lower, double upper) {
            return lower + (upper - lower) * uniform_dist(*rng);
        }

        // Ratios at least one are accepted without consuming another random number.
        bool accept(double ratio) {
            return ratio >= 1.0 || uniform_dist(*rng) < ratio;
        }

        /** @brief Boltzmann factor for a change already integrated over imaginary time. */
        static double boltzmann_weight(double energy_diff) {
            if (energy_diff == 0.) return 1.;
            return std::exp(-energy_diff);
        }

        /** @brief Estimate tau_int in samples; warn when it exceeds 10% of the series length. */
        static double calculate_autocorrelation_time_with_warning(
            const std::vector<double>& obs_real,
            const std::string& observable_name,
            bool has_hysteresis_context = false,
            size_t hysteresis_point = 0,
            double h = 0.,
            double lmbda = 0.
        );

        /**
         * @brief Recompute the coupled diagonal energy integral from event histories.
         *
         * Returns -h * edge_integral - mu * star_integral in the x-basis, or
         * -lmbda * edge_integral - J * plaquette_integral in the z-basis.
         * Zero-coupling terms are skipped; the running energy and caches are unchanged.
         */
        static double total_integrated_pot_energy(Lattice& lat, double h, double mu, double J, double lmbda);

        /**
         * @brief Rebuild diagonal caches, rotate the time origin, and recompute the energy.
         *
         * Call before each stage that changes couplings, since updates may skip caches
         * whose coupling is zero. Periodic calls also limit accumulated rounding error.
         * The global time rotation preserves full-period integrals.
         */
        static void reinitialize_potential_energy(
            Lattice& lat, double& integrated_pot_energy,
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Edge contribution for reversing one spin on an ordered time interval.
         *
         * @return (coupled change, bare edge-integral change). The coupled value is
         *         -h times the bare change in the x-basis, or -lmbda times it in z.
         * @param total_cache Use -2 * cached integral only for a full-period flip.
         * @pre imag_time_spin_flip < imag_time_next_spin_flip.
         * @note A zero active edge coupling returns two zeros without reading the cache.
         */
        static std::tuple<double, double> 
        integrated_pot_energy_diff_single_spin_flip_edge(
            Lattice& lat, double h, double mu, double J, double lmbda, 
            const Lattice::Edge& edg, double imag_time_spin_flip, double imag_time_next_spin_flip, 
            bool total_cache
        );

        /**
         * @brief Diagonal tuple contribution for reversing one edge's spin.
         *
         * @return (coupled change, affected indices, aligned bare integral changes).
         *         Indices are endpoint stars in x, or adjacent plaquettes in z.
         *         Only the scalar is multiplied by -mu (x) or -J (z).
         * @param total_cache Use -2 * cached integrals only for a full-period flip.
         * @pre imag_time_spin_flip < imag_time_next_spin_flip.
         * @note A zero active tuple coupling returns zero and empty vectors.
         */
        static std::tuple<double, SmallIndexVector, SmallEnergyVector> 
        integrated_pot_energy_diff_single_spin_flip_tuple(
            Lattice& lat, double h, double mu, double J, double lmbda, 
            const Lattice::Edge& edg, double imag_time_spin_flip, double imag_time_next_spin_flip, 
            bool total_cache
        );

        /**
         * @brief Edge contributions for reversing every spin of an update tuple.
         *
         * @return (coupled sum, borrowed tuple_edges view, aligned bare changes).
         *         Only the sum is multiplied by -h (x) or -lmbda (z).
         * @param total_cache Use -2 * cached integrals only for a full-period flip.
         * @param interval_has_no_inner_flips Use the constant-spin fast path; the caller
         *        must ensure no event lies strictly inside the interval on any edge.
         * @param known_flip_indices Optional ranks of known_flip_time in each edge's
         *        full history, aligned with tuple_edges; empty requests fresh searches.
         * @param known_flip_time Event time bordering the interval for the cached ranks.
         * @pre imag_time_spin_flip < imag_time_next_spin_flip.
         * @note A zero active coupling returns a zero sum and one zero per edge.
         */
        static std::tuple<double, std::span<const Lattice::Edge>, SmallEnergyVector> 
        integrated_pot_energy_diff_tuple_flip_edge(
            Lattice& lat, double h, double mu, double J, double lmbda, 
            int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
            double imag_time_spin_flip, double imag_time_next_spin_flip, 
            bool total_cache,
            bool interval_has_no_inner_flips = false,
            std::span<const int> known_flip_indices = {}, double known_flip_time = 0.
        );

        /**
         * @brief Edge contributions for a tuple event paired with one single event per edge.
         *
         * @return (coupled sum, borrowed tuple_edges view, aligned bare changes).
         *         Only the sum is multiplied by -h (x) or -lmbda (z).
         * @param imag_time_spin_flips Single-event times aligned with tuple_edges.
         * @param create_vector Creation/deletion flags; unused when computing the integral.
         * @param tuple_destroy Tuple creation/deletion flag; unused in the integral.
         * @pre tau_left < tau_right and imag_time_spin_flips.size() == tuple_edges.size().
         * @note Equal-time event pairs cancel by parity. A zero active coupling returns
         *       a zero sum and one zero per edge.
         */
        static std::tuple<double, std::span<const Lattice::Edge>, SmallEnergyVector> 
        integrated_pot_energy_diff_combination_flip_edge(
            Lattice& lat, double h, double mu, double J, double lmbda, 
            int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
            double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
            double tau_left, double tau_right, 
            const SmallBoolVector& create_vector, bool tuple_destroy
        );

        /**
         * @brief Diagonal tuple contributions for a combination update.
         *
         * @return (coupled sum, sorted affected indices, aligned bare changes).
         *         Affected tuples are stars in x or plaquettes in z; only the scalar
         *         is multiplied by -mu (x) or -J (z).
         * @param tuple_index Update plaquette index in x, or star center in z.
         * @param imag_time_spin_flips Single-event times aligned with tuple_edges.
         * @param create_vector Creation/deletion flags; unused in the integral.
         * @param tuple_destroy Tuple creation/deletion flag; unused in the integral.
         * @pre tau_left < tau_right and imag_time_spin_flips.size() == tuple_edges.size().
         * @note A zero active tuple coupling returns zero and empty vectors.
         */
        static std::tuple<double, SmallIndexVector, SmallEnergyVector> 
        integrated_pot_energy_diff_combination_flip_tuple(
            Lattice& lat, double h, double mu, double J, double lmbda, 
            int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
            double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
            double tau_left, double tau_right, 
            const SmallBoolVector& create_vector, bool tuple_destroy
        );

        /**
         * @brief Commit the event-history changes of an accepted combination update.
         *
         * @param tuple_index Plaquette index in x, or star center in z.
         * @param tuple_edges Edges of that tuple, in the order of the per-edge arrays.
         * @param imag_time_tuple_flip Tuple event to insert or remove.
         * @param imag_time_spin_flips Single-event times aligned with tuple_edges.
         * @param create_vector True inserts a single event; false removes an existing one.
         * @param tuple_destroy True removes the tuple event; false inserts it.
         * @pre Both per-edge arrays have tuple_edges.size() entries.
         * @note The caller updates energy caches separately. h and mu are unused.
         */
        static void combination_flip(
            Lattice& lat, double h, double mu, 
            int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
            double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
            const SmallBoolVector& create_vector, bool tuple_destroy
        );

        /**
         * @brief Propose inserting or removing two single-spin events on one edge.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_double_single_spin_flip(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Propose moving one single-spin event within its neighboring-event window.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_single_spin_flip_move(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Propose reversing one edge's spin over the full imaginary-time period.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_global_single_spin_flip(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );


        /**
         * @brief Propose reversing every spin of one tuple over the full time period.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_global_tuple_flip(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Propose inserting or removing two events on one update tuple.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_double_tuple_flip(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Propose moving a tuple event within the common window of its edges.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_single_tuple_flip_move(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Propose creating/removing a tuple event and one single event per edge.
         * @see metropolis_step() for shared argument and acceptance-diagnostic contracts.
         */
        void metropolis_step_spin_tuple_combination(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        /**
         * @brief Choose one of the seven proposal types with equal probability.
         *
         * @param lat Lattice whose histories and active caches are changed on acceptance.
         * @param integrated_pot_energy Running coupled integral; changed only on acceptance.
         * @param acc_ratio Receives the raw Metropolis ratio, possibly greater than one,
         *        or zero for an abandoned proposal. Acceptance uses min(1, ratio);
         *        this diagnostic is not an accepted/rejected flag.
         * @param beta Imaginary-time period.
         * @param h Electric-field coupling.
         * @param mu Star coupling.
         * @param J Plaquette coupling.
         * @param lmbda Gauge-field coupling.
         */
        void metropolis_step(
            Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
            double h, double mu, double J, double lmbda
        );

        template<typename T>
        requires std::integral<T> || std::floating_point<T>
        constexpr T modulo(T a, T b) {
            if constexpr (std::integral<T>) {
                T result = a % b;
                return result >= 0 ? result : result + b;
            } else {
                T result = std::fmod(a, b);
                return result >= 0 ? result : result + b;
            }
        }

        template<std::floating_point T>
        constexpr bool almost_equal(
            T a, T b, 
            T rel_tol = std::numeric_limits<T>::epsilon(), 
            T abs_tol = std::numeric_limits<T>::min()
        ) {
            T diff = std::abs(a - b);

            if (diff <= abs_tol) return true;

            T largest = std::max(std::abs(a), std::abs(b));
            return diff <= largest * rel_tol;
        }
};

template<char Basis>
requires ValidBasis<Basis>
std::vector<std::function<std::variant< std::complex<double>, double>(Lattice&, double, double, double, double)>> 
ExtendedToricCodeQMC<Basis>::get_obs_func_vec(const std::vector<std::string>& observables) {
    std::vector<std::function<std::variant< std::complex<double>, double>(Lattice&, double, double, double, double)>> result;
    for (const auto& obs_name : observables) {
        bool obs_found = false;
        for (const auto& obs : obs_vec) {
            if (obs.obs_name == obs_name) {
                result.emplace_back( obs.obs_func );
                obs_found = true;
                break;
            } 
        }
        if (!obs_found) {
            throw std::invalid_argument(std::string("The observable \"") + obs_name + std::string("\" does not exist!"));
        }
    }
    return result;
}

template<char Basis>
requires ValidBasis<Basis>
std::vector<std::string> 
ExtendedToricCodeQMC<Basis>::get_obs_type_vec(const std::vector<std::string>& observables) {
    std::vector<std::string> result;
    for (const auto& obs_name : observables) {
        bool obs_found = false;
        for (const auto& obs : obs_vec) {
            if (obs.obs_name == obs_name) {
                result.emplace_back( obs.obs_type );
                obs_found = true;
                break;
            } 
        }
        if (!obs_found) {
            throw std::invalid_argument(std::string("The observable \"") + obs_name + std::string("\" does not exist!"));
        }
    }
    return result;
}

template<char Basis>
requires ValidBasis<Basis>
double ExtendedToricCodeQMC<Basis>::total_integrated_pot_energy(
    Lattice& lat, double h, double mu, double J, double lmbda
) {
    if constexpr (Basis == 'x') {
        if (h == 0.) {
            if (mu == 0.) return 0.;
            return -lat.total_integrated_star_energy() * mu;
        }
        if (mu == 0.) return -lat.total_integrated_edge_energy() * h;
        return -lat.total_integrated_edge_energy() * h - lat.total_integrated_star_energy() * mu;
    } else if constexpr (Basis == 'z') {
        if (lmbda == 0.) {
            if (J == 0.) return 0.;
            return -lat.total_integrated_plaquette_energy() * J;
        }
        if (J == 0.) return -lat.total_integrated_edge_energy() * lmbda;
        return -lat.total_integrated_edge_energy() * lmbda - lat.total_integrated_plaquette_energy() * J;
    }
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::reinitialize_potential_energy(
    Lattice& lat, double& integrated_pot_energy,
    double h, double mu, double J, double lmbda
) {
    lat.init_potential_energy();
    lat.rotate_imag_time();
    integrated_pot_energy = total_integrated_pot_energy(lat, h, mu, J, lmbda);
}

template<char Basis>
requires ValidBasis<Basis>
std::tuple<double, double> 
ExtendedToricCodeQMC<Basis>::integrated_pot_energy_diff_single_spin_flip_edge(
    Lattice& lat, double h, double mu, double J, double lmbda, 
    const Lattice::Edge& edg, double imag_time_spin_flip, double imag_time_next_spin_flip, 
    bool total_cache
) {
    // Disabled terms need neither an acceptance contribution nor a cache update.
    // Parameter-changing workflows rebuild the cache before re-enabling them.
    if constexpr (Basis == 'x') {
        if (h == 0.) return {0., 0.};
    } else if constexpr (Basis == 'z') {
        if (lmbda == 0.) return {0., 0.};
    }

    double delta_energy_single = 0.;
    double bare_energy_single = 0.;
    if (total_cache) {
        bare_energy_single = -2*lat.get_potential_edge_energy(edg);
    } else {
        bare_energy_single = lat.integrated_edge_energy_diff(edg, imag_time_spin_flip, imag_time_next_spin_flip);
    }

    if constexpr (Basis == 'x') {
        delta_energy_single = -h * bare_energy_single;
    } else if constexpr (Basis == 'z') {
        delta_energy_single = -lmbda * bare_energy_single;
    }

    return {delta_energy_single, bare_energy_single};
}

template<char Basis>
requires ValidBasis<Basis>
std::tuple<double, typename ExtendedToricCodeQMC<Basis>::SmallIndexVector, typename ExtendedToricCodeQMC<Basis>::SmallEnergyVector> 
ExtendedToricCodeQMC<Basis>::integrated_pot_energy_diff_single_spin_flip_tuple(
    Lattice& lat, double h, double mu, double J, double lmbda, 
    const Lattice::Edge& edg, double imag_time_spin_flip, double imag_time_next_spin_flip, 
    bool total_cache
) {
    double delta_energy_tuple = 0.;
    // Disabled terms need neither an acceptance contribution nor a cache update.
    // Parameter-changing workflows rebuild the cache before re-enabling them.
    if constexpr (Basis == 'x') {
        if (mu == 0.) return {0., SmallIndexVector{}, SmallEnergyVector{}};
        auto [bare_energy, star_centers, bare_star_potential_energy_diffs]
        = lat.integrated_star_energy_diff(edg, imag_time_spin_flip, imag_time_next_spin_flip, total_cache);
        delta_energy_tuple = -mu * bare_energy; // The bare change already includes the factor -2.
        return {delta_energy_tuple, std::move(star_centers), std::move(bare_star_potential_energy_diffs)};
    } else if constexpr (Basis == 'z') {
        if (J == 0.) return {0., SmallIndexVector{}, SmallEnergyVector{}};
        auto [bare_energy, plaquette_indices, bare_plaquette_potential_energy_diffs]
        = lat.integrated_plaquette_energy_diff(edg, imag_time_spin_flip, imag_time_next_spin_flip, total_cache);
        delta_energy_tuple = -J * bare_energy; // The bare change already includes the factor -2.
        return {delta_energy_tuple, std::move(plaquette_indices), std::move(bare_plaquette_potential_energy_diffs)};
    } 
    return {0., SmallIndexVector{}, SmallEnergyVector{}};
}

template<char Basis>
requires ValidBasis<Basis>
std::tuple<double, std::span<const Lattice::Edge>, typename ExtendedToricCodeQMC<Basis>::SmallEnergyVector> 
ExtendedToricCodeQMC<Basis>::integrated_pot_energy_diff_tuple_flip_edge(
    Lattice& lat, double h, double mu, double J, double lmbda, 
    int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
    double imag_time_spin_flip, double imag_time_next_spin_flip, 
    bool total_cache,
    bool interval_has_no_inner_flips,
    std::span<const int> known_flip_indices, double known_flip_time
) {
    UNUSED(tuple_index);
    UNUSED(mu);
    UNUSED(J);
    // Keep a correctly-sized zero delta vector for accepted tuple updates.
    if constexpr (Basis == 'x') {
        if (h == 0.) return {0., tuple_edges, SmallEnergyVector(tuple_edges.size(), 0.)};
    } else if constexpr (Basis == 'z') {
        if (lmbda == 0.) return {0., tuple_edges, SmallEnergyVector(tuple_edges.size(), 0.)};
    }

    double delta_energy_single = 0.;
    SmallEnergyVector bare_energy_single_vector;
    bare_energy_single_vector.reserve(tuple_edges.size());
    double bare_energy_edg = 0.;
    for (size_t edge_index = 0; edge_index < tuple_edges.size(); ++edge_index) {
        const Lattice::Edge& edg = tuple_edges[edge_index];
        if (total_cache) { 
            bare_energy_edg = -2*lat.get_potential_edge_energy(edg);
        } else if (interval_has_no_inner_flips) {
            bare_energy_edg = lat.integrated_edge_energy_diff_no_inner_flips(
                edg, imag_time_spin_flip, imag_time_next_spin_flip,
                known_flip_indices.empty() ? -1 : known_flip_indices[edge_index], known_flip_time
            );
        } else {
            bare_energy_edg = lat.integrated_edge_energy_diff(edg, imag_time_spin_flip, imag_time_next_spin_flip);
        }
        bare_energy_single_vector.emplace_back(bare_energy_edg);

        if constexpr (Basis == 'x') {
            delta_energy_single += -h * bare_energy_edg;
        } else if constexpr (Basis == 'z') {
            delta_energy_single += -lmbda * bare_energy_edg;
        }
    }
    return {delta_energy_single, tuple_edges, std::move(bare_energy_single_vector)};
}

template<char Basis>
requires ValidBasis<Basis>
std::tuple<double, std::span<const Lattice::Edge>, typename ExtendedToricCodeQMC<Basis>::SmallEnergyVector> 
ExtendedToricCodeQMC<Basis>::integrated_pot_energy_diff_combination_flip_edge(
    Lattice& lat, double h, double mu, double J, double lmbda, 
    int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
    double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
    double tau_left, double tau_right, 
    const SmallBoolVector& create_vector, bool tuple_destroy
) {
    // Keep a correctly-sized zero delta vector for accepted combination updates.
    if constexpr (Basis == 'x') {
        if (h == 0.) return {0., tuple_edges, SmallEnergyVector(tuple_edges.size(), 0.)};
    } else if constexpr (Basis == 'z') {
        if (lmbda == 0.) return {0., tuple_edges, SmallEnergyVector(tuple_edges.size(), 0.)};
    }

    double energy_single_diff = 0.; 
    SmallEnergyVector bare_energy_single_vector;
    bare_energy_single_vector.reserve(tuple_edges.size());
    for (size_t i = 0; i < tuple_edges.size(); ++i) {
        const auto edg = tuple_edges[i];

        double bare_energy_edg = lat.integrated_edge_energy_diff_combination(
            edg,
            tau_left,
            tau_right,
            imag_time_tuple_flip,
            imag_time_spin_flips[i]
        );
        bare_energy_single_vector.emplace_back(bare_energy_edg);
        if constexpr (Basis == 'x') {
            energy_single_diff += -h * bare_energy_edg;
        }  else if constexpr (Basis == 'z') {
            energy_single_diff += -lmbda * bare_energy_edg;
        }
    }

    return {energy_single_diff, tuple_edges, std::move(bare_energy_single_vector)};
}

template<char Basis>
requires ValidBasis<Basis>
std::tuple<double, typename ExtendedToricCodeQMC<Basis>::SmallIndexVector, typename ExtendedToricCodeQMC<Basis>::SmallEnergyVector> 
ExtendedToricCodeQMC<Basis>::integrated_pot_energy_diff_combination_flip_tuple(
    Lattice& lat, double h, double mu, double J, double lmbda, 
    int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
    double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
    double tau_left, double tau_right, 
    const SmallBoolVector& create_vector, bool tuple_destroy
) {
    double energy_tuple_diff = 0.; 
    // Disabled terms need neither an acceptance contribution nor a cache update.
    if constexpr (Basis == 'x') {
        if (mu == 0.) return {0., SmallIndexVector{}, SmallEnergyVector{}};
        auto [bare_energy, star_centers, bare_star_potential_energy_diffs]
        = lat.integrated_star_energy_diff_combination(
            tuple_index, tau_left, tau_right, 
            std::span<const double>(imag_time_spin_flips.data(), imag_time_spin_flips.size()), 
            imag_time_tuple_flip
        );
        energy_tuple_diff += -mu * bare_energy;
        return {energy_tuple_diff, std::move(star_centers), std::move(bare_star_potential_energy_diffs)};
    }  else if constexpr (Basis == 'z') {
        if (J == 0.) return {0., SmallIndexVector{}, SmallEnergyVector{}};
        auto [bare_energy, plaquette_indices, bare_plaquette_potential_energy_diffs] =
        lat.integrated_plaquette_energy_diff_combination(
            tuple_index, tau_left, tau_right, 
            std::span<const double>(imag_time_spin_flips.data(), imag_time_spin_flips.size()), 
            imag_time_tuple_flip
        );
        energy_tuple_diff += -J * bare_energy;
        return {energy_tuple_diff, std::move(plaquette_indices), std::move(bare_plaquette_potential_energy_diffs)};
    }
    return {0., SmallIndexVector{}, SmallEnergyVector{}};
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::combination_flip(
    Lattice& lat, double h, double mu, 
    int tuple_index, std::span<const Lattice::Edge> tuple_edges, 
    double imag_time_tuple_flip, const SmallEnergyVector& imag_time_spin_flips, 
    const SmallBoolVector& create_vector, bool tuple_destroy
) {
    UNUSED(h);
    UNUSED(mu);
    if (tuple_destroy) {
        for (size_t i = 0; i < tuple_edges.size(); ++i) {
            const auto edg = tuple_edges[i];

            if (create_vector[i]) {
                lat.insert_single_spin_flip(edg, imag_time_spin_flips[i]);
            } else {
                lat.delete_single_spin_flip(edg, lat.get_spin_flip_index(edg, imag_time_spin_flips[i]));
            }
        }
        lat.delete_tuple_flip(tuple_index, tuple_edges, imag_time_tuple_flip);
    } else {
        for (size_t i = 0; i < tuple_edges.size(); ++i) {
            const auto edg = tuple_edges[i];

            if (create_vector[i]) {
                lat.insert_single_spin_flip(edg, imag_time_spin_flips[i]);
            } else {
                lat.delete_single_spin_flip(edg, lat.get_spin_flip_index(edg, imag_time_spin_flips[i]));
            }
        }
        lat.insert_tuple_flip(tuple_index, tuple_edges, imag_time_tuple_flip);
    } 
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_double_single_spin_flip(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_double_single_spin_flip!";
#endif

    if (Basis == 'x' && (lmbda == 0)) {
        acc_ratio = 0.; 
        return;
    }

    if (Basis == 'z' && (h == 0)) {
        acc_ratio = 0.; 
        return;
    }

    const double rnd_create_destroy = uniform_dist(*rng);

#ifndef NDEBUG
    const auto& [rand_edge, source, target] = lat.get_random_edge();
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- Randomly chose edge between vertices {} and {}.", source, target);
#else
    const auto rand_edge = lat.get_random_edge_descriptor();
#endif
    
    std::span<const double> single_spin_flips = lat.get_single_spin_flips(rand_edge);
    int single_spin_flip_count = single_spin_flips.size();

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- Found {} single spin flips.", single_spin_flip_count);
#endif

    if (rnd_create_destroy < 0.5) { // destroys a pair of single spin flips
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_single_spin_flip --- Trying to DESTROY single spin flip pair.";
#endif
        if (single_spin_flip_count > 1) [[likely]] {
            const int random_spin_flip_index = random_index(single_spin_flip_count - 1);

            const double imag_time_spin_flip = single_spin_flips[random_spin_flip_index];

            const double imag_time_next_spin_flip = single_spin_flips[random_spin_flip_index+1];

            double imag_time_next_next_spin_flip = beta;

            if (random_spin_flip_index != single_spin_flip_count - 2) [[likely]] {
                imag_time_next_next_spin_flip = single_spin_flips[random_spin_flip_index+2];
            }

            double imag_time_prev_spin_flip = 0.;

            if (random_spin_flip_index != 0) [[likely]] {
                imag_time_prev_spin_flip = single_spin_flips[random_spin_flip_index-1];
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- Imaginary time of random single spin flip: {}; previous: {}; next: {}; next next: {}", 
                imag_time_spin_flip, imag_time_prev_spin_flip, imag_time_next_spin_flip, imag_time_next_next_spin_flip);
#endif            
            const auto& [integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge] 
            = integrated_pot_energy_diff_single_spin_flip_edge(
                lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, imag_time_next_spin_flip, false
            );
            const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
            = integrated_pot_energy_diff_single_spin_flip_tuple(
                lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, imag_time_next_spin_flip, false
            );
            const double integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;
            if constexpr (Basis == 'x') {
                acc_ratio = 2./(lmbda * lmbda * (imag_time_next_next_spin_flip - imag_time_prev_spin_flip) 
                * (imag_time_next_next_spin_flip - imag_time_prev_spin_flip)) * boltzmann_weight(integrated_pot_energy_diff);
            } else if constexpr (Basis == 'z') {
                acc_ratio = 2./(h * h * (imag_time_next_next_spin_flip - imag_time_prev_spin_flip) 
                * (imag_time_next_next_spin_flip - imag_time_prev_spin_flip)) * boltzmann_weight(integrated_pot_energy_diff);
            }
            

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- acceptance ratio: {}", acc_ratio);
#endif   
            
            if (accept(acc_ratio)) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_single_spin_flip --- ACCEPTED.";
#endif   
                lat.delete_double_single_spin_flip(rand_edge, imag_time_spin_flip, imag_time_next_spin_flip);
                integrated_pot_energy += integrated_pot_energy_diff;
                lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                    int tuple_index = pot_energy_tuple_indices[i];
                    double potential_tuple_energy_diff = pot_energy_diffs[i];
                    if constexpr (Basis == 'x') {
                        lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                    } else if constexpr (Basis == 'z') {
                        lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                    }
                }
            } 
        } else [[unlikely]] {  
            acc_ratio = 0.;
        }
    } else { // create a pair of single spin flips
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_single_spin_flip --- Trying to CREATE single spin flip pair.";
#endif
        double tau_left, tau_right;
        int random_spin_flip_index = 0;

        if (single_spin_flip_count == 0) {
            tau_left = 0.;
            tau_right = beta;
        } else {
            random_spin_flip_index = random_index(single_spin_flip_count + 1);

            if (random_spin_flip_index == 0) [[unlikely]] {
                tau_left = 0.; // not really necessary
                tau_right = single_spin_flips[0];
            } else if (random_spin_flip_index == single_spin_flip_count) [[unlikely]] {
                tau_left = single_spin_flips[random_spin_flip_index-1];
                tau_right = beta;
            } else [[likely]] {
                tau_left = single_spin_flips[random_spin_flip_index-1];
                tau_right = single_spin_flips[random_spin_flip_index];
            }
        }

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- tau_left: {}; tau_right: {}", tau_left, tau_right);
#endif

        double tau_1 = 0., tau_2 = 0.;

        tau_1 = uniform_real(tau_left, tau_right);
        tau_2 = uniform_real(tau_left, tau_right);

        if (tau_2 < tau_1) {
            std::swap(tau_1, tau_2);
        }

        if (std::abs(tau_2 - tau_1) < PRECISION 
        || std::abs(tau_1 - tau_left) < PRECISION 
        || std::abs(tau_2 - tau_right) < PRECISION
        ) [[unlikely]] {
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_single_spin_flip --- Random numbers very close to each other.";
#endif
            acc_ratio = 0.; 
            return;
        } else [[likely]] {
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- tau_1: {}; tau_2: {}", tau_1, tau_2);
#endif
            const auto& [integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge] 
            = integrated_pot_energy_diff_single_spin_flip_edge(lat, h, mu, J, lmbda, rand_edge, tau_1, tau_2, false);
            const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
            = integrated_pot_energy_diff_single_spin_flip_tuple(lat, h, mu, J, lmbda, rand_edge, tau_1, tau_2, false);
            const double integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;
            if constexpr (Basis == 'x') {
                acc_ratio = ((tau_right - tau_left)*(tau_right - tau_left))/2. 
                * lmbda * lmbda * boltzmann_weight(integrated_pot_energy_diff);
            } else if constexpr (Basis == 'z') {
                acc_ratio = ((tau_right - tau_left)*(tau_right - tau_left))/2. 
                * h * h * boltzmann_weight(integrated_pot_energy_diff);
            }
            

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_single_spin_flip --- acceptance ratio: {}", acc_ratio);
#endif   

            if (accept(acc_ratio)) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_single_spin_flip --- ACCEPTED.";
#endif  
                lat.insert_double_single_spin_flip(rand_edge, tau_1, tau_2);
                integrated_pot_energy += integrated_pot_energy_diff;
                lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                    int tuple_index = pot_energy_tuple_indices[i];
                    double potential_tuple_energy_diff = pot_energy_diffs[i];
                    if constexpr (Basis == 'x') {
                        lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                    } else if constexpr (Basis == 'z') {
                        lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                    }
                }
            } 
        }  
    }
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_single_spin_flip_move(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_single_spin_flip_move!";
#endif

    if (Basis == 'x' && (lmbda == 0)) {
        acc_ratio = 0.; 
        return;
    }

    if (Basis == 'z' && (h == 0)) {
        acc_ratio = 0.; 
        return;
    }

#ifndef NDEBUG
    const auto& [rand_edge, source, target] = lat.get_random_edge();
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- Randomly chose edge between vertices {} and {}.", source, target);    
#else
    const auto rand_edge = lat.get_random_edge_descriptor();
#endif

    std::span<const double> single_spin_flips = lat.get_single_spin_flips(rand_edge);
    int single_spin_flip_count = single_spin_flips.size();

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- Found {} single spin flips.", single_spin_flip_count);  
#endif    

    if (single_spin_flip_count != 0) [[likely]] {
        const int random_spin_flip_index = random_index(single_spin_flip_count);

        const double imag_time_spin_flip = single_spin_flips[random_spin_flip_index];

        const double imag_time_next_spin_flip = lat.flip_next_imag_time(rand_edge, imag_time_spin_flip);

        const double imag_time_prev_spin_flip = lat.flip_prev_imag_time(rand_edge, imag_time_spin_flip);

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- Imaginary time of random single spin flip: {}; previous: {}; next: {}", imag_time_spin_flip, imag_time_prev_spin_flip, imag_time_next_spin_flip);
#endif  

        const int random_spin_flip_index_lat = lat.get_spin_flip_index(rand_edge, imag_time_spin_flip);

        if (imag_time_prev_spin_flip < imag_time_spin_flip && imag_time_spin_flip < imag_time_next_spin_flip) {
            const double new_imag_time = uniform_real(imag_time_prev_spin_flip, imag_time_next_spin_flip);

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- New imaginary time: {}", new_imag_time);
#endif  

            if (std::abs(new_imag_time - imag_time_spin_flip) < PRECISION) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            } else {
                double integrated_pot_energy_diff = 0.;
                double integrated_pot_energy_diff_edge = 0.;
                double bare_pot_energy_diff_edge = 0.;
                SmallIndexVector pot_energy_tuple_indices;
                SmallEnergyVector pot_energy_diffs;
                if (new_imag_time > imag_time_spin_flip) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices_tmp, pot_energy_diffs_tmp]
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    pot_energy_tuple_indices = std::move(pot_energy_tuple_indices_tmp);
                    pot_energy_diffs = std::move(pot_energy_diffs_tmp);
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;
                } else {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices_tmp, pot_energy_diffs_tmp]
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    pot_energy_tuple_indices = std::move(pot_energy_tuple_indices_tmp);
                    pot_energy_diffs = std::move(pot_energy_diffs_tmp);
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;
                }
                acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                if (accept(acc_ratio)) {
#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED.";
#endif              
                    lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, true);
                    integrated_pot_energy += integrated_pot_energy_diff;
                    lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                    for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                        int tuple_index = pot_energy_tuple_indices[i];
                        double potential_tuple_energy_diff = pot_energy_diffs[i];
                        if constexpr (Basis == 'x') {
                            lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                        } else if constexpr (Basis == 'z') {
                            lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                        }
                    }
                } 

            }
        } else { // Potentially cross beta
            const double new_imag_time = modulo(uniform_real(imag_time_prev_spin_flip, beta + imag_time_next_spin_flip), beta);

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- New imaginary time: {}", new_imag_time);
#endif  
            
            if (std::abs(new_imag_time - imag_time_spin_flip) < PRECISION || std::abs(new_imag_time - imag_time_prev_spin_flip) < PRECISION || std::abs(new_imag_time - imag_time_next_spin_flip) < PRECISION) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            } else {
                double integrated_pot_energy_diff = 0.;
                double integrated_pot_energy_diff_edge = 0.;
                double bare_pot_energy_diff_edge = 0.;

                if (random_spin_flip_index_lat == 0 && new_imag_time < imag_time_spin_flip) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED.";
#endif 
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, true);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                } else if (random_spin_flip_index_lat > 0 && imag_time_spin_flip < new_imag_time) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED.";
#endif 
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, true);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                } else if (random_spin_flip_index_lat == 0 && imag_time_spin_flip < new_imag_time && new_imag_time < imag_time_next_spin_flip) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, new_imag_time, false
                    );
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED.";
#endif 
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, true);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                } else if (random_spin_flip_index_lat > 0 
                && imag_time_spin_flip > new_imag_time && new_imag_time > imag_time_prev_spin_flip
                ) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, imag_time_spin_flip, false
                    );
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED.";
#endif 
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, true);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                } else if (random_spin_flip_index_lat == 0 
                && imag_time_spin_flip < new_imag_time && new_imag_time > imag_time_next_spin_flip
                ) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, beta, false
                    );
                    auto [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, new_imag_time, beta, false
                    );
                    if (imag_time_spin_flip > 0) {
                        auto tmp = integrated_pot_energy_diff_single_spin_flip_edge(
                            lat, h, mu, J, lmbda, rand_edge, 0., imag_time_spin_flip, false
                        );
                        integrated_pot_energy_diff_edge += std::get<0>(tmp);
                        bare_pot_energy_diff_edge += std::get<1>(tmp);
                        const auto& [integrated_pot_energy_diff_tuple_2, pot_energy_tuple_indices_2, pot_energy_diffs_2] 
                        = integrated_pot_energy_diff_single_spin_flip_tuple(
                            lat, h, mu, J, lmbda, rand_edge, 0., imag_time_spin_flip, false
                        );
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            pot_energy_diffs[i] += pot_energy_diffs_2[i];
                        }
                        integrated_pot_energy_diff_tuple += integrated_pot_energy_diff_tuple_2;
                    }
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED. Flipping spin.";
#endif 
                        lat.flip_spin(rand_edge);
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, false);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                } else if (random_spin_flip_index_lat > 0 
                && imag_time_spin_flip > new_imag_time && new_imag_time < imag_time_prev_spin_flip
                ) {
                    std::tie(integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge) 
                    = integrated_pot_energy_diff_single_spin_flip_edge(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, beta, false
                    );
                    auto [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
                    = integrated_pot_energy_diff_single_spin_flip_tuple(
                        lat, h, mu, J, lmbda, rand_edge, imag_time_spin_flip, beta, false
                    );
                    if (new_imag_time > 0) {
                        auto tmp = integrated_pot_energy_diff_single_spin_flip_edge(
                            lat, h, mu, J, lmbda, rand_edge, 0., new_imag_time, false
                        );
                        integrated_pot_energy_diff_edge += std::get<0>(tmp);
                        bare_pot_energy_diff_edge += std::get<1>(tmp);
                        const auto& [integrated_pot_energy_diff_tuple_2, pot_energy_tuple_indices_2, pot_energy_diffs_2] 
                        = integrated_pot_energy_diff_single_spin_flip_tuple(
                            lat, h, mu, J, lmbda, rand_edge, 0., new_imag_time, false
                        );
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            pot_energy_diffs[i] += pot_energy_diffs_2[i];
                        }
                        integrated_pot_energy_diff_tuple += integrated_pot_energy_diff_tuple_2;
                    }
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_spin_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_spin_flip_move --- ACCEPTED. Flipping spin.";
#endif 
                        lat.flip_spin(rand_edge);
                        lat.move_spin_flip(rand_edge, random_spin_flip_index_lat, new_imag_time, false);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
                        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
                            int tuple_index = pot_energy_tuple_indices[i];
                            double potential_tuple_energy_diff = pot_energy_diffs[i];
                            if constexpr (Basis == 'x') {
                                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
                            } else if constexpr (Basis == 'z') {
                                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
                            }
                        }
                    } 
                }
            }
        }
    }
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_global_single_spin_flip(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_global_single_spin_flip!";
#endif

#ifndef NDEBUG
    const auto& [rand_edge, source, target] = lat.get_random_edge();
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_global_single_spin_flip --- Randomly chose edge between vertices {} and {}.", source, target);
#else
    const auto rand_edge = lat.get_random_edge_descriptor();
#endif

    const auto& [integrated_pot_energy_diff_edge, bare_pot_energy_diff_edge] 
    = integrated_pot_energy_diff_single_spin_flip_edge(
        lat, h, mu, J, lmbda, rand_edge, 0., beta, true
    );
    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_diffs] 
    = integrated_pot_energy_diff_single_spin_flip_tuple(
        lat, h, mu, J, lmbda, rand_edge, 0., beta, true
    );
    const double integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_global_single_spin_flip --- acceptance ratio: {}", acc_ratio);
#endif  

    if (accept(acc_ratio)) {
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_global_single_spin_flip --- ACCEPTED.";
#endif 
        lat.flip_spin(rand_edge);
        integrated_pot_energy += integrated_pot_energy_diff;
        lat.add_potential_edge_energy(rand_edge, bare_pot_energy_diff_edge);
        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
            int tuple_index = pot_energy_tuple_indices[i];
            double potential_tuple_energy_diff = pot_energy_diffs[i];
            if constexpr (Basis == 'x') {
                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
            } else if constexpr (Basis == 'z') {
                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
            }
        } 
    } 
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_global_tuple_flip(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_global_tuple_flip!";
#endif

    int random_tuple = -1;
    std::span<const Lattice::Edge> tuple_edges;

    if constexpr (Basis == 'x') {
        random_tuple = lat.get_random_plaquette_index();
        tuple_edges = lat.get_plaquette_edges(random_tuple);
    } else if constexpr (Basis == 'z') {
        random_tuple = lat.get_random_vertex();
        tuple_edges = lat.get_star_edges(random_tuple);
    }
    if (tuple_edges.empty()) return;

#ifndef NDEBUG
    if constexpr (Basis == 'x') {
        const auto tuple_vertex_pairs = lat.get_plaquette_vertex_pairs(random_tuple);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_global_tuple_flip --- Randomly chose plaquette with index {} and edges {}", random_tuple, tuple_vertex_pairs);
    } else if constexpr (Basis == 'z') {
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_global_tuple_flip --- Randomly chose star with center {}", random_tuple);
    }
#endif

    const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
    = integrated_pot_energy_diff_tuple_flip_edge(lat, h, mu, J, lmbda, random_tuple, tuple_edges, 0, beta, true);

    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_global_single_spin_flip --- acceptance ratio: {}", acc_ratio);
#endif  

    if (accept(acc_ratio)) {
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_global_tuple_flip --- ACCEPTED.";
#endif 
        integrated_pot_energy += integrated_pot_energy_diff;
        for (size_t i = 0; i < tuple_edges.size(); ++i) {
            auto edg = tuple_edges[i];
            lat.flip_spin(edg);
            lat.add_potential_edge_energy(edg, pot_energy_diffs[i]);
        }
    } 

}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_double_tuple_flip(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_double_tuple_flip!";
#endif

    if (Basis == 'x' && (J == 0)) {
        acc_ratio = 0.; 
        return;
    }

    if (Basis == 'z' && (mu == 0)) {
        acc_ratio = 0.; 
        return;
    }

    const double rnd_create_destroy = uniform_dist(*rng);
    
    int random_tuple = -1;
    std::span<const Lattice::Edge> tuple_edges;

    if constexpr (Basis == 'x') {
        random_tuple = lat.get_random_plaquette_index();
        tuple_edges = lat.get_plaquette_edges(random_tuple);
    } else if constexpr (Basis == 'z') {
        random_tuple = lat.get_random_vertex();
        tuple_edges = lat.get_star_edges(random_tuple);
    }
    if (tuple_edges.empty()) return;

#ifndef NDEBUG
    if constexpr (Basis == 'x') {
        const auto tuple_vertex_pairs = lat.get_plaquette_vertex_pairs(random_tuple);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- Randomly chose plaquette with index {} and edges {}", random_tuple, tuple_vertex_pairs);
    } else if constexpr (Basis == 'z') {
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- Randomly chose star with center ", random_tuple);
    }
#endif

    std::span<const double> tuple_flips = lat.get_tuple_spin_flips(random_tuple);
    int tuple_flip_count = tuple_flips.size();

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- Found {} tuple flips.", tuple_flip_count);
#endif  

    if (rnd_create_destroy < 0.5) { // destroy tuples
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_tuple_flip --- Trying to DESTROY tuple flip pair.";
#endif
        if (tuple_flip_count > 1) [[likely]] {
            const int random_tuple_flip_index = random_index(tuple_flip_count - 1);

            const double imag_time_tuple_flip = tuple_flips[random_tuple_flip_index];

            const double imag_time_next_tuple_flip = tuple_flips[random_tuple_flip_index+1];

            double imag_time_next_next_tuple_flip = beta;

            if (random_tuple_flip_index != tuple_flip_count - 2) [[likely]] {
                imag_time_next_next_tuple_flip = tuple_flips[random_tuple_flip_index+2];
            }

            double imag_time_prev_tuple_flip = 0.;

            if (random_tuple_flip_index != 0) [[likely]] {
                imag_time_prev_tuple_flip = tuple_flips[random_tuple_flip_index-1];
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- Imaginary time of random tuple flip: {}; previous: {}; next: {}; next next: {}", imag_time_tuple_flip, imag_time_prev_tuple_flip, imag_time_next_tuple_flip, imag_time_next_next_tuple_flip);
#endif         

            const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
            = integrated_pot_energy_diff_tuple_flip_edge(
                lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                imag_time_tuple_flip, imag_time_next_tuple_flip, false
            );

            if constexpr (Basis == 'x') {
                acc_ratio = 2./(J * J * (imag_time_next_next_tuple_flip - imag_time_prev_tuple_flip) 
                * (imag_time_next_next_tuple_flip - imag_time_prev_tuple_flip)) 
                * boltzmann_weight(integrated_pot_energy_diff);
            } else if constexpr (Basis == 'z') {
                acc_ratio = 2./(mu * mu * (imag_time_next_next_tuple_flip - imag_time_prev_tuple_flip) 
                * (imag_time_next_next_tuple_flip - imag_time_prev_tuple_flip)) 
                * boltzmann_weight(integrated_pot_energy_diff);
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- acceptance ratio: ", acc_ratio);
#endif  
            
            if (accept(acc_ratio)) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_tuple_flip --- ACCEPTED.";
#endif  
                lat.delete_double_tuple_flip(
                    random_tuple, tuple_edges, imag_time_tuple_flip, imag_time_next_tuple_flip
                );
                integrated_pot_energy += integrated_pot_energy_diff;
                for (size_t i = 0; i < tuple_edges.size(); ++i) {
                    lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                }
            }
        } else [[unlikely]] {
            acc_ratio = 0.;
        }
    } else { // create tuples
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_tuple_flip --- Trying to CREATE tuple flip pair.";
#endif
        double tau_left, tau_right;
        int random_tuple_flip_index = 0;

        if (tuple_flip_count > 0) [[likely]] {
            random_tuple_flip_index = random_index(tuple_flip_count + 1);

            if (random_tuple_flip_index == 0) {
                tau_left = 0.; //not really necessary
                tau_right = tuple_flips[0];
            } else if (random_tuple_flip_index == tuple_flip_count) {
                tau_left = tuple_flips[random_tuple_flip_index-1];
                tau_right = beta;
            } else {
                tau_left = tuple_flips[random_tuple_flip_index-1];
                tau_right = tuple_flips[random_tuple_flip_index];
            }
        } else [[unlikely]] {
            tau_left = 0.;
            tau_right = beta;
        }

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- tau_left: {}; tau_right: {}", tau_left, tau_right);
#endif

        double tau_1 = 0., tau_2 = 0.;

        tau_1 = uniform_real(tau_left, tau_right);
        tau_2 = uniform_real(tau_left, tau_right);
        if (tau_2 < tau_1) {
            std::swap(tau_1, tau_2);
        }

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- tau_1: {}; tau_2: {}", tau_1, tau_2);
#endif

        if (std::abs(tau_2 - tau_1) < PRECISION 
        || std::abs(tau_1 - tau_left) < PRECISION 
        || std::abs(tau_2 - tau_right) < PRECISION
        ) [[unlikely]] {
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_tuple_flip --- Random numbers very close to each other.";
#endif
            acc_ratio = 0.; 
            return;
        } else {
            const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
            = integrated_pot_energy_diff_tuple_flip_edge(
                lat, h, mu, J, lmbda, random_tuple, tuple_edges, tau_1, tau_2, false
            );

            if constexpr (Basis == 'x') {
                acc_ratio = ((tau_right - tau_left)*(tau_right - tau_left))/2. 
                * J * J * boltzmann_weight(integrated_pot_energy_diff);
            } else if constexpr (Basis == 'z') {
                acc_ratio = ((tau_right - tau_left)*(tau_right - tau_left))/2. 
                * mu * mu * boltzmann_weight(integrated_pot_energy_diff);
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_double_tuple_flip --- acceptance ratio: {}", acc_ratio);
#endif 

            if (accept(acc_ratio)) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_double_tuple_flip --- ACCEPTED.";
#endif  
                lat.insert_double_tuple_flip(random_tuple, tuple_edges, tau_1, tau_2);
                integrated_pot_energy += integrated_pot_energy_diff;
                for (size_t i = 0; i < tuple_edges.size(); ++i) {
                    lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                }
            } 
        }
    }
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_single_tuple_flip_move(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_single_tuple_flip_move!";
#endif

    if (Basis == 'x' && (J == 0)) {
        acc_ratio = 0.; 
        return;
    }

    if (Basis == 'z' && (mu == 0)) {
        acc_ratio = 0.; 
        return;
    }

    int random_tuple = -1;
    std::span<const Lattice::Edge> tuple_edges;

    if constexpr (Basis == 'x') {
        random_tuple = lat.get_random_plaquette_index();
        tuple_edges = lat.get_plaquette_edges(random_tuple);
    } else if constexpr (Basis == 'z') {
        random_tuple = lat.get_random_vertex();
        tuple_edges = lat.get_star_edges(random_tuple);
    }
    if (tuple_edges.empty()) [[unlikely]] return;

#ifndef NDEBUG
    if constexpr (Basis == 'x') {
        const auto tuple_vertex_pairs = lat.get_plaquette_vertex_pairs(random_tuple);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- Randomly chose plaquette with index {} and edges {}", random_tuple, tuple_vertex_pairs);
    } else if constexpr (Basis == 'z') {
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- Randomly chose star with center {}", random_tuple);
    }
#endif

    std::span<const double> tuple_flips = lat.get_tuple_spin_flips(random_tuple);
    int tuple_flip_count = tuple_flips.size();

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- Found {} tuple flips.", tuple_flip_count);
#endif  

    if (tuple_flip_count > 0) [[likely]] {
        const int random_tuple_flip_index = random_index(tuple_flip_count);

        const double imag_time_tuple_flip = tuple_flips[random_tuple_flip_index];

        SmallIndexVector edge_flip_indices;
        const auto [tau_left, tau_right] =
            lat.tuple_flip_window(tuple_edges, imag_time_tuple_flip, &edge_flip_indices);

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("ropolis_step_single_tuple_flip_move --- Imaginary time of random tuple flip: {}; tau_left: {}; tau_right: {}", imag_time_tuple_flip, tau_left, tau_right);
#endif  

        if (tau_left < imag_time_tuple_flip && imag_time_tuple_flip < tau_right) { // not cross beta
            const double new_imag_time = uniform_real(tau_left, tau_right);

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- New imaginary time: {}", new_imag_time);
#endif 

            if (std::abs(new_imag_time - imag_time_tuple_flip) < PRECISION) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            } else {
                double integrated_pot_energy_diff = 0.;
                std::span<const Lattice::Edge> pot_energy_edges;
                SmallEnergyVector pot_energy_diffs;
                if (new_imag_time > imag_time_tuple_flip) {
                    const auto& [integrated_pot_energy_diff_edge, pot_energy_edges_tmp, pot_energy_diffs_tmp]
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        imag_time_tuple_flip, new_imag_time, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    pot_energy_edges = pot_energy_edges_tmp;
                    pot_energy_diffs = std::move(pot_energy_diffs_tmp);
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge;
                } else {
                    const auto& [integrated_pot_energy_diff_edge, pot_energy_edges_tmp, pot_energy_diffs_tmp]
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        new_imag_time, imag_time_tuple_flip, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    pot_energy_edges = pot_energy_edges_tmp;
                    pot_energy_diffs = std::move(pot_energy_diffs_tmp);
                    integrated_pot_energy_diff = integrated_pot_energy_diff_edge;
                }
                acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                if (accept(acc_ratio)) {
#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED.";
#endif  
                    lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, true,
                                        edge_flip_indices, random_tuple_flip_index);
                    integrated_pot_energy += integrated_pot_energy_diff;
                    for (size_t i = 0; i < tuple_edges.size(); ++i) {
                        lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                    }
                } 

            }
        } else { // Potentially cross beta
            const double new_imag_time = modulo(uniform_real(tau_left, beta + tau_right), beta);

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- New imaginary time: {}", new_imag_time);
#endif 


            if (std::abs(new_imag_time - imag_time_tuple_flip) < PRECISION 
            || std::abs(new_imag_time - tau_left) < PRECISION 
            || std::abs(new_imag_time - tau_right) < PRECISION
            ) {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            } else if (tau_left > imag_time_tuple_flip) { // Potentially over beta left
                if (new_imag_time < imag_time_tuple_flip) { // Normal move
                    const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        new_imag_time, imag_time_tuple_flip, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED.";
#endif  
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, true,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    } 
                } else if (new_imag_time > imag_time_tuple_flip && new_imag_time < tau_right) { // Normal move
                    const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        imag_time_tuple_flip, new_imag_time, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED.";
#endif  
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, true,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    } 
                } else if (new_imag_time > imag_time_tuple_flip && new_imag_time > tau_left) { // Move over beta
                    auto [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        new_imag_time, beta, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    if (imag_time_tuple_flip != 0) {
                        const auto& [integrated_pot_energy_diff_2, pot_energy_edges_2, pot_energy_diffs_2] 
                        = integrated_pot_energy_diff_tuple_flip_edge(
                            lat, h, mu, J, lmbda, random_tuple, tuple_edges,
                             0., imag_time_tuple_flip, false, true, edge_flip_indices, imag_time_tuple_flip
                        );
                        for (size_t i = 0; i < pot_energy_edges.size(); ++i) {
                            pot_energy_diffs[i] += pot_energy_diffs_2[i];
                        }
                        integrated_pot_energy_diff += integrated_pot_energy_diff_2;
                    }

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED. Flipping spins.";
#endif  
                        for (const auto& p_edg : tuple_edges) {
                            lat.flip_spin(p_edg);
                        }
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, false,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    }
                }
            } else if (tau_right < imag_time_tuple_flip) { // Potentially over beta right
                if (new_imag_time > imag_time_tuple_flip) { // Normal move
                    const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        imag_time_tuple_flip, new_imag_time, false, true, edge_flip_indices, imag_time_tuple_flip
                    );

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED.";
#endif  
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, true,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    } 
                } else if (new_imag_time < imag_time_tuple_flip && new_imag_time > tau_left) { // Normal move
                    const auto& [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        new_imag_time, imag_time_tuple_flip, false, true, edge_flip_indices, imag_time_tuple_flip
                    );

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED.";
#endif  
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, true,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    } 
                } else if (new_imag_time < imag_time_tuple_flip && new_imag_time < tau_right) { // Move over beta
                    auto [integrated_pot_energy_diff, pot_energy_edges, pot_energy_diffs] 
                    = integrated_pot_energy_diff_tuple_flip_edge(
                        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                        imag_time_tuple_flip, beta, false, true, edge_flip_indices, imag_time_tuple_flip
                    );
                    if (new_imag_time > 0) {
                        const auto& [integrated_pot_energy_diff_2, pot_energy_edges_2, pot_energy_diffs_2] 
                        = integrated_pot_energy_diff_tuple_flip_edge(
                            lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
                            0., new_imag_time, false, true, edge_flip_indices, imag_time_tuple_flip
                        );
                        for (size_t i = 0; i < pot_energy_edges.size(); ++i) {
                            pot_energy_diffs[i] += pot_energy_diffs_2[i];
                        }
                        integrated_pot_energy_diff += integrated_pot_energy_diff_2;
                    }

                    acc_ratio = boltzmann_weight(integrated_pot_energy_diff);

#ifndef NDEBUG
                    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- acceptance ratio: {}", acc_ratio);
#endif  

                    if (accept(acc_ratio)) {
#ifndef NDEBUG
                        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_single_tuple_flip_move --- ACCEPTED. Flipping spins.";
#endif 
                        for (const auto& p_edg : tuple_edges) {
                            lat.flip_spin(p_edg);
                        }
                        lat.move_tuple_flip(random_tuple, tuple_edges, imag_time_tuple_flip, new_imag_time, false,
                                            edge_flip_indices, random_tuple_flip_index);
                        integrated_pot_energy += integrated_pot_energy_diff;
                        for (size_t i = 0; i < tuple_edges.size(); ++i) {
                            lat.add_potential_edge_energy(tuple_edges[i], pot_energy_diffs[i]);
                        }
                    } 
                }
            }
        }
    }
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step_spin_tuple_combination(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << "";
    BOOST_LOG_TRIVIAL(debug) << "# Welcome to MC update metropolis_step_spin_tuple_combination!";
#endif

    if (Basis == 'x' && (J == 0 || lmbda == 0)) {
        acc_ratio = 0.; 
        return;
    }

    if (Basis == 'z' && (mu == 0 || h == 0)) {
        acc_ratio = 0.; 
        return;
    }
    
    double rnd_create_destroy = uniform_dist(*rng);
    
    int random_tuple = -1;
    std::span<const Lattice::Edge> tuple_edges;

    if constexpr (Basis == 'x') {
        random_tuple = lat.get_random_plaquette_index();
        tuple_edges = lat.get_plaquette_edges(random_tuple);
    } else if constexpr (Basis == 'z') {
        random_tuple = lat.get_random_vertex();
        tuple_edges = lat.get_star_edges(random_tuple);
    }
    if (tuple_edges.empty()) return;

#ifndef NDEBUG
    if constexpr (Basis == 'x') {
        const auto tuple_vertex_pairs = lat.get_plaquette_vertex_pairs(random_tuple);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- Randomly chose plaquette with index {} and edges {}", random_tuple, tuple_vertex_pairs);
    } else if constexpr (Basis == 'z') {
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_single_tuple_flip_move --- Randomly chose star with center {}", random_tuple);
    }
#endif

    std::span<const double> tuple_flips = lat.get_tuple_spin_flips(random_tuple);
    int tuple_flip_count = tuple_flips.size();

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- Found {} tuple flips.", tuple_flip_count);
#endif  

    bool tuple_destroy = true;
    int random_tuple_flip_index = -1;
    double tau_new = -1., r = 0., tau_left = 0., tau_right = 0.;

    if (rnd_create_destroy < 0.5) { // Destroys tuple
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_spin_tuple_combination --- Trying to DESTROY tuple flip.";
#endif
        tuple_destroy = true;
        if (tuple_flip_count > 0) [[likely]] {
            random_tuple_flip_index = random_index(tuple_flip_count);

            double imag_time_tuple_flip = tuple_flips[random_tuple_flip_index];
            tau_new = imag_time_tuple_flip;

            double imag_time_next_tuple_flip = beta;

            if (random_tuple_flip_index < tuple_flip_count-1) {
                imag_time_next_tuple_flip = tuple_flips[(random_tuple_flip_index+1)];
            }

            tau_right = imag_time_next_tuple_flip;

            double imag_time_prev_tuple_flip = 0.;

            if (random_tuple_flip_index != 0) {
                imag_time_prev_tuple_flip = tuple_flips[(random_tuple_flip_index-1)];
            }
            tau_left = imag_time_prev_tuple_flip;

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- Imaginary time of random tuple flip: {}; tau_left: {}; tau_right: {}", imag_time_tuple_flip, tau_left, tau_right);
#endif  
            if constexpr (Basis == 'x') {
                r = 1. / ((imag_time_next_tuple_flip-imag_time_prev_tuple_flip) * J);
            } else if constexpr (Basis == 'z') {
                r = 1. / ((imag_time_next_tuple_flip-imag_time_prev_tuple_flip) * mu);
            }
        } else [[unlikely]] {
            acc_ratio = 0.; 
            return;
        }
    } else { // Create tuple
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_spin_tuple_combination --- Trying to CREATE tuple flip.";
#endif
        tuple_destroy = false;

        random_tuple_flip_index = 0;

        if (tuple_flip_count > 0) [[likely]] {
            random_tuple_flip_index = random_index(tuple_flip_count + 1);

            if (random_tuple_flip_index == 0) {
                tau_left = 0.; // not really necessary
                tau_right = tuple_flips[0];
            } else if (random_tuple_flip_index == tuple_flip_count) {
                tau_left = tuple_flips[random_tuple_flip_index-1];
                tau_right = beta;
            } else {
                tau_left = tuple_flips[random_tuple_flip_index-1];
                tau_right = tuple_flips[random_tuple_flip_index];
            }
        } else [[unlikely]] {
            tau_left = 0.;
            tau_right = beta;
        }

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- tau_left: {}; tau_right: {}", tau_left, tau_right);
#endif  

        if (tau_left < tau_right) {
            tau_new = uniform_real(tau_left, tau_right);

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- Proposed tuple flip imaginary time: {}", tau_new);
#endif 

#ifndef NDEBUG
            if (lat.check_spin_flips_present_tuple(tuple_edges, tau_new)) {
                acc_ratio = 0.; 
                return;
            }
#endif

            if (std::abs(tau_new - tau_left) < PRECISION || std::abs(tau_new - tau_right) < PRECISION) [[unlikely]] {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_spin_tuple_combination --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            }
            
            if constexpr (Basis == 'x') {
                r = (tau_right-tau_left) * J;
            } else if constexpr (Basis == 'z') {
                r = (tau_right-tau_left) * mu;
            }
        } else {
            acc_ratio = 0.; 
            return;
        }
    }

    SmallEnergyVector r_b;
    SmallBoolVector create_vector;
    SmallEnergyVector flip_times;
    r_b.reserve(tuple_edges.size());
    create_vector.reserve(tuple_edges.size());
    flip_times.reserve(tuple_edges.size());

    // Bound the time interval affected by the tuple and all per-edge events.
    double tau_right_potential_energy = tau_right;
    double tau_left_potential_energy = tau_left;

    for (size_t i = 0; i < tuple_edges.size(); ++i) {
        rnd_create_destroy = uniform_dist(*rng);
        int tau_spin_flip_index = -1;
        double tau_spin_flip = -1.;
        
        const Lattice::Edge edg = tuple_edges[i];
        
        std::span<const double> single_spin_flips = lat.get_single_spin_flips(edg);
        int single_spin_flip_count = single_spin_flips.size();
        const auto next_flip_it = detail::time_upper_bound(
            single_spin_flips.begin(), single_spin_flips.end(), tau_new
        );
        const int imag_time_next_flip_index = static_cast<int>(
            next_flip_it - single_spin_flips.begin()
        );
        const int imag_time_prev_flip_index = imag_time_next_flip_index - 1;
        const double imag_time_prev_flip = imag_time_prev_flip_index >= 0
            ? single_spin_flips[imag_time_prev_flip_index]
            : 0.;
        const double imag_time_next_flip = next_flip_it != single_spin_flips.end()
            ? *next_flip_it
            : beta;

#ifndef NDEBUG
        int spin = lat.get_spin(edg);
        const auto& [source_v, target_v] = lat.vertices_of_edge(edg);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - Spin: {}. Found {} single spin flips.", source_v, target_v, spin, single_spin_flip_count);
#endif  

        if (rnd_create_destroy < 0.5) { // destroy
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - Trying to DESTROY single spin flip.", source_v, target_v);
#endif  

            if (single_spin_flip_count == 0) [[unlikely]] {
                acc_ratio = 0.; 
                return;
            }

            create_vector.emplace_back(false);
            rnd_create_destroy = uniform_dist(*rng);

            if (rnd_create_destroy < 0.5) {
                tau_spin_flip = imag_time_prev_flip;
                tau_spin_flip_index = imag_time_prev_flip_index;
            } else {
                tau_spin_flip = imag_time_next_flip;
                tau_spin_flip_index = imag_time_next_flip_index;
            }

            if (tau_spin_flip_index == -1 || tau_spin_flip_index == single_spin_flip_count) [[unlikely]] {
                acc_ratio = 0.; 
                return;
            }
            
            double tau_prev_of_chosen = 0.;
            if (tau_spin_flip_index > 0) {
                tau_prev_of_chosen = single_spin_flips[tau_spin_flip_index-1];
            }

            double tau_next_of_chosen = beta;
            if (tau_spin_flip_index < single_spin_flip_count-1) {
                tau_next_of_chosen = single_spin_flips[tau_spin_flip_index+1];
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - Propose to remove single spin flip at {}", source_v, target_v, tau_spin_flip);            
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - previous: {} - next: {}", source_v, target_v, tau_prev_of_chosen, tau_next_of_chosen);
#endif 

            if constexpr (Basis == 'x') {
                r_b.emplace_back(2./((tau_next_of_chosen - tau_prev_of_chosen) * lmbda));
            } else if constexpr (Basis == 'z') {
                r_b.emplace_back(2./((tau_next_of_chosen - tau_prev_of_chosen) * h));
            }

            flip_times.emplace_back(tau_spin_flip);

            if (tau_next_of_chosen > tau_right_potential_energy) {
                tau_right_potential_energy = tau_next_of_chosen;
            }

            if (tau_prev_of_chosen < tau_left_potential_energy) {
                tau_left_potential_energy = tau_prev_of_chosen;
            }
        } else { // create
#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - Trying to CREATE single spin flip.", source_v, target_v);
#endif  
            
            create_vector.emplace_back(true);

            tau_spin_flip = uniform_real(imag_time_prev_flip, imag_time_next_flip);

            if (std::abs(tau_spin_flip - imag_time_prev_flip) < PRECISION 
            || std::abs(tau_spin_flip - imag_time_next_flip) < PRECISION
            ) [[unlikely]] {
#ifndef NDEBUG
                BOOST_LOG_TRIVIAL(debug) << "metropolis_step_spin_tuple_combination --- Random numbers very close to each other.";
#endif
                acc_ratio = 0.; 
                return;
            }

            if (std::abs(tau_spin_flip - tau_new) < PRECISION) [[unlikely]] {
                acc_ratio = 0.; 
                return;
            }

#ifndef NDEBUG
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - Propose to create single spin flip at {}", source_v, target_v, tau_spin_flip);
            BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- edge between vertices {} and {} - previous: {} - next: {}", source_v, target_v, imag_time_prev_flip, imag_time_next_flip);
#endif 
            
            if constexpr (Basis == 'x') {
                r_b.emplace_back((imag_time_next_flip - imag_time_prev_flip) * lmbda / 2.);
            } else if constexpr (Basis == 'z') {
                r_b.emplace_back((imag_time_next_flip - imag_time_prev_flip) * h / 2.);
            }

            flip_times.emplace_back(tau_spin_flip);

            if (imag_time_next_flip > tau_right_potential_energy) {
                tau_right_potential_energy = imag_time_next_flip;
            }

            if (imag_time_prev_flip < tau_left_potential_energy) {
                tau_left_potential_energy = imag_time_prev_flip;
            }     
        }
    }

    const auto& [integrated_pot_energy_diff_edge, pot_energy_edges, pot_energy_edges_diffs] 
    = integrated_pot_energy_diff_combination_flip_edge(
        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
        tau_new, flip_times, tau_left_potential_energy, tau_right_potential_energy, 
        create_vector, tuple_destroy
    );
    const auto& [integrated_pot_energy_diff_tuple, pot_energy_tuple_indices, pot_energy_tuple_diffs] 
    = integrated_pot_energy_diff_combination_flip_tuple(
        lat, h, mu, J, lmbda, random_tuple, tuple_edges, 
        tau_new, flip_times, tau_left_potential_energy, tau_right_potential_energy, 
        create_vector, tuple_destroy
    );
    const double integrated_pot_energy_diff = integrated_pot_energy_diff_edge + integrated_pot_energy_diff_tuple;

#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- integrated_pot_energy_diff_edge: {}", integrated_pot_energy_diff_edge);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- integrated_pot_energy_diff_tuple: {}", integrated_pot_energy_diff_tuple);
        BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- integrated_pot_energy_diff: {}", integrated_pot_energy_diff);
#endif 

    acc_ratio = boltzmann_weight(integrated_pot_energy_diff) * r * std::accumulate(r_b.begin(), r_b.end(), 1., std::multiplies<double>());

#ifndef NDEBUG
    BOOST_LOG_TRIVIAL(debug) << std::format("metropolis_step_spin_tuple_combination --- acceptance ratio: {}", acc_ratio);
#endif  

    if (accept(acc_ratio)) {
#ifndef NDEBUG
        BOOST_LOG_TRIVIAL(debug) << "metropolis_step_spin_tuple_combination --- ACCEPTED.";
#endif 
        combination_flip(lat, h, mu, random_tuple, tuple_edges, tau_new, flip_times, create_vector, tuple_destroy); 
        integrated_pot_energy += integrated_pot_energy_diff;
        for (size_t i = 0; i < tuple_edges.size(); ++i) {
            lat.add_potential_edge_energy(pot_energy_edges[i], pot_energy_edges_diffs[i]);
        }
        for (size_t i = 0; i < pot_energy_tuple_indices.size(); ++i) {
            int tuple_index = pot_energy_tuple_indices[i];
            double potential_tuple_energy_diff = pot_energy_tuple_diffs[i];
            if constexpr (Basis == 'x') {
                lat.add_potential_star_energy(tuple_index, potential_tuple_energy_diff);
            } else if constexpr (Basis == 'z') {
                lat.add_potential_plaquette_energy(tuple_index, potential_tuple_energy_diff);
            }
        } 
    } 
}

template<char Basis>
requires ValidBasis<Basis>
void ExtendedToricCodeQMC<Basis>::metropolis_step(
    Lattice& lat, double& integrated_pot_energy, double& acc_ratio, double beta, 
    double h, double mu, double J, double lmbda
) {
    const int rnd = random_index(7);
    if (rnd < 1) {
        metropolis_step_double_single_spin_flip(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else if (rnd < 2) {
        metropolis_step_single_spin_flip_move(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else if (rnd < 3) {
        metropolis_step_global_single_spin_flip(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else if (rnd < 4) {
        metropolis_step_global_tuple_flip(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else if (rnd < 5) {
        metropolis_step_double_tuple_flip(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else if (rnd < 6) {
        metropolis_step_single_tuple_flip_move(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    } else {
        metropolis_step_spin_tuple_combination(lat, integrated_pot_energy, acc_ratio, beta, h, mu, J, lmbda);
    }  

#ifndef NDEBUG
    double integrated_pot_energy_check = total_integrated_pot_energy(
        lat, h, mu, J, lmbda
    );

    BOOST_LOG_TRIVIAL(debug) << std::format("Integrated potential energy --- Cached: {}, Actual: {}", integrated_pot_energy, integrated_pot_energy_check);
#endif  
}

template<char Basis>
requires ValidBasis<Basis>
double ExtendedToricCodeQMC<Basis>::calculate_autocorrelation_time_with_warning(
    const std::vector<double>& obs_real,
    const std::string& observable_name,
    bool has_hysteresis_context,
    size_t hysteresis_point,
    double h,
    double lmbda
) {
    const double autocorrelation_time = paratoric::statistics::get_autocorrelation_time(
        paratoric::statistics::get_autocorrelation_function(obs_real)
    );
    const size_t sample_count = obs_real.size();
    const double warning_threshold =
        AUTOCORRELATION_WARNING_SAMPLE_FRACTION * static_cast<double>(sample_count);

    if (sample_count > 0 && autocorrelation_time > warning_threshold) {
        const std::string hysteresis_context = has_hysteresis_context
            ? std::format(" at hysteresis point {} (h={}, lmbda={})", hysteresis_point, h, lmbda)
            : "";

        BOOST_LOG_TRIVIAL(warning)
            << "Autocorrelation time for observable \""
            << observable_name << "\"" << hysteresis_context
            << " is " << autocorrelation_time
            << ", greater than " << AUTOCORRELATION_WARNING_SAMPLE_FRACTION
            << " * number of samples (" << warning_threshold
            << " for " << sample_count << " samples).";
    }

    return autocorrelation_time;
}

template<char Basis>
requires ValidBasis<Basis>
Result ExtendedToricCodeQMC<Basis>::get_thermalization(
    const Config& config
) {

    detail::validate_thermalization_config(config);

    if constexpr (Basis != 'x' && Basis != 'z') {
        throw std::invalid_argument("Basis must be either \"x\" or \"z\".");
    }

    if (Basis != config.lat_spec.basis) {
        throw std::invalid_argument("Template parameter basis and config.lat_spec basis must match.");
    }
    
    if constexpr (Basis == 'x') {
        if (config.param_spec.J < 0) {
            throw std::invalid_argument("J must be non-negative in the x-basis.");
        } else if (config.param_spec.lmbda < 0) {
            throw std::invalid_argument("lmbda must be non-negative in the x-basis.");
        }
    } else if constexpr (Basis == 'z') {
        if (config.param_spec.mu < 0) {
            throw std::invalid_argument("mu must be non-negative in the z-basis.");
        } else if (config.param_spec.h < 0) {
            throw std::invalid_argument("h must be non-negative in the z-basis.");
        }
    }

    if (config.sim_spec.seed != 0) rng->set_seed(config.sim_spec.seed);

    const auto thermalization_count = static_cast<size_t>(std::max(config.sim_spec.N_thermalization, 0));
    std::vector<double> acc_ratio_vector;
    acc_ratio_vector.reserve(thermalization_count);
    // Series: [observable][sample], in the requested observable order.
    std::vector<std::vector<std::variant< std::complex<double>, double>>> observable_vector;
    observable_vector.reserve(config.sim_spec.observables.size());
    for (const auto& obs_func : config.sim_spec.observables) {
        UNUSED(obs_func);
        observable_vector.emplace_back();
        observable_vector.back().reserve(thermalization_count);
    } 

    auto obs_func_vec = get_obs_func_vec(config.sim_spec.observables);

    auto lat = Lattice(config.lat_spec, rng);
    
    double integrated_pot_energy = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );
    double acc_ratio = 1.;

    int metropolis_step_count = 0;
    int reset_potential_energy_count = static_cast<int>(lat.get_edge_count()*10000);

    for (int i = 0; i < config.sim_spec.N_thermalization; ++i) {
        ++metropolis_step_count;
        metropolis_step(
            lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, 
            config.param_spec.h, config.param_spec.mu, 
            config.param_spec.J, config.param_spec.lmbda
        );

        if (metropolis_step_count == reset_potential_energy_count) [[unlikely]] {
            metropolis_step_count = 0;
            // Rebuild caches periodically to limit accumulated rounding error.
            reinitialize_potential_energy(
                lat, integrated_pot_energy, config.param_spec.h, config.param_spec.mu,
                config.param_spec.J, config.param_spec.lmbda
            );
        }

        for (size_t k = 0; k < config.sim_spec.observables.size(); k++) {
            observable_vector[k].emplace_back(
                obs_func_vec[k](lat, config.param_spec.h, config.param_spec.lmbda, config.param_spec.mu, config.param_spec.J)
            );
        }
        acc_ratio_vector.emplace_back(acc_ratio);

        if (config.out_spec.save_snapshots && i%10000 == 0) [[unlikely]] {
            lat.update_spin_string();
        }
    }

    double integrated_pot_energy_check = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );

    if (!almost_equal(integrated_pot_energy, integrated_pot_energy_check, 1e-5, 1e-13)) {
        throw std::runtime_error(std::format("Integrated potential energy mismatch. {} does not match {}.", integrated_pot_energy, integrated_pot_energy_check));
    }

    if (config.out_spec.save_snapshots) {
        lat.write_graph("snapshots", config.out_spec.path_out);
    }

    return Result{
        .series=std::move(observable_vector), 
        .acc_ratio=std::move(acc_ratio_vector)
    };                                     
}

template<char Basis>
requires ValidBasis<Basis>
Result ExtendedToricCodeQMC<Basis>::get_sample(
    const Config& config
) { 

    detail::validate_sample_config(config);

    if constexpr (Basis != 'x' && Basis != 'z') {
        throw std::invalid_argument("Basis must be either \"x\" or \"z\".");
    }

    if (Basis != config.lat_spec.basis) {
        throw std::invalid_argument("Template parameter basis and config.lat_spec basis must match.");
    }

    if constexpr (Basis == 'x') {
        if (config.param_spec.J < 0) {
            throw std::invalid_argument("J must be non-negative in the x-basis.");
        } else if (config.param_spec.lmbda < 0) {
            throw std::invalid_argument("lmbda must be non-negative in the x-basis.");
        }
    } else if constexpr (Basis == 'z') {
        if (config.param_spec.mu < 0) {
            throw std::invalid_argument("mu must be non-negative in the z-basis.");
        } else if (config.param_spec.h < 0) {
            throw std::invalid_argument("h must be non-negative in the z-basis.");
        }
    }

    if (config.sim_spec.seed != 0) rng->set_seed(config.sim_spec.seed);

    auto obs_func_vec = get_obs_func_vec(config.sim_spec.observables);
    auto obs_type_vec = get_obs_type_vec(config.sim_spec.observables);
    
    // Series: [observable][sample], in the requested observable order.
    const auto sample_count = static_cast<size_t>(std::max(config.sim_spec.N_samples, 0));
    std::vector<std::vector< std::variant< std::complex<double>, double> >> observable_vector;
    observable_vector.reserve(config.sim_spec.observables.size());
    std::vector<double> observable_mean_vector(config.sim_spec.observables.size(), 0.), 
                        observable_std_vector(config.sim_spec.observables.size(), 0.), 
                        binder_mean_vector(config.sim_spec.observables.size(), 0.), 
                        binder_std_vector(config.sim_spec.observables.size(), 0.), 
                        observable_autocorrelation_time_vector(config.sim_spec.observables.size(), 0.);
    for (const auto& obs_func : config.sim_spec.observables) {
        UNUSED(obs_func);
        observable_vector.emplace_back();
        observable_vector.back().reserve(sample_count);
    } 

    auto lat = Lattice(config.lat_spec, rng);
    
    double integrated_pot_energy = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );
    double acc_ratio = 1.;

    if (config.sim_spec.custom_therm) {
        integrated_pot_energy = total_integrated_pot_energy(
            lat, config.param_spec.h_therm, config.param_spec.mu,
            config.param_spec.J, config.param_spec.lmbda_therm
        );

        // Prepare the state at the custom fields before ramping to the target.
        for (int i = 0; i < config.sim_spec.N_thermalization; ++i) {
            metropolis_step(
                lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, 
                config.param_spec.h_therm, config.param_spec.mu, config.param_spec.J, 
                config.param_spec.lmbda_therm
            );
        }

        // Ramp h first, then lmbda; rebuild caches at every change of fields.
        double h_end = config.param_spec.h_therm;
        if (std::abs(config.param_spec.h - config.param_spec.h_therm) > PRECISION) {
            for (int i = 9; i > -1; --i) {
                h_end = config.param_spec.h + (config.param_spec.h_therm-config.param_spec.h) * i / 10.;
                reinitialize_potential_energy(
                    lat, integrated_pot_energy, h_end, config.param_spec.mu,
                    config.param_spec.J, config.param_spec.lmbda_therm
                );
                for (int j = 0; j < config.sim_spec.N_thermalization / 10.; ++j) {
                    metropolis_step(
                        lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, h_end, 
                        config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda_therm
                    );
                }
            }
        }
        double lmbda_end = config.param_spec.lmbda_therm;
        if (std::abs(config.param_spec.lmbda - config.param_spec.lmbda_therm) > PRECISION) {
            for (int i = 9; i > -1; --i) {
                lmbda_end = config.param_spec.lmbda + (config.param_spec.lmbda_therm-config.param_spec.lmbda) * i / 10.;
                reinitialize_potential_energy(
                    lat, integrated_pot_energy, h_end, config.param_spec.mu,
                    config.param_spec.J, lmbda_end
                );
                for (int j = 0; j < config.sim_spec.N_thermalization / 10.; ++j) {
                    metropolis_step(
                        lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, h_end, 
                        config.param_spec.mu, config.param_spec.J, lmbda_end
                    );
                }
            }
        }
    } else {
        // Thermalization 
        for (int i = 0; i < config.sim_spec.N_thermalization; ++i) {
            metropolis_step(
                lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, 
                config.param_spec.h, config.param_spec.mu, config.param_spec.J, 
                config.param_spec.lmbda
            );
        }
    }

    double integrated_pot_energy_check = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );

    if (!almost_equal(integrated_pot_energy, integrated_pot_energy_check, 1e-5, 1e-13)) {
        throw std::runtime_error(std::format("Integrated potential energy mismatch. {} does not match {}.", integrated_pot_energy, integrated_pot_energy_check));
    }

    int metropolis_step_count = 0;
    int reset_potential_energy_count = static_cast<int>(lat.get_edge_count()*100000);

    for (int i = 0; i < config.sim_spec.N_samples; ++i) {
        for (int j = 0; j < config.sim_spec.N_between_samples; ++j) {
            ++metropolis_step_count;
            metropolis_step(
                lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, config.param_spec.h, 
                config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
            );
            if (metropolis_step_count == reset_potential_energy_count) [[unlikely]] {
                metropolis_step_count = 0;
                // Rebuild caches periodically to limit accumulated rounding error.
                reinitialize_potential_energy(
                    lat, integrated_pot_energy, config.param_spec.h, config.param_spec.mu,
                    config.param_spec.J, config.param_spec.lmbda
                );
            }
        }

        for (size_t k = 0; k < config.sim_spec.observables.size(); k++) {
            observable_vector[k].emplace_back(
                obs_func_vec[k](lat, config.param_spec.h, config.param_spec.lmbda, config.param_spec.mu, config.param_spec.J)
            );
        }

        if (config.out_spec.save_snapshots) {
            lat.update_spin_string();
        }
    }

    if (config.out_spec.save_snapshots) {
        lat.write_graph("snapshots", config.out_spec.path_out);
    }

    for (size_t k = 0; k < config.sim_spec.observables.size(); k++) {
        if (obs_type_vec[k] == "real") {
            std::vector<double> obs_real;
            const auto& series = observable_vector[k];
            obs_real.reserve(series.size());
            for (auto const& x : series) {
                obs_real.emplace_back(std::get<double>(x));
            }

            const auto& [observable_mean, observable_std, binder_mean, binder_std] 
            = paratoric::statistics::get_bootstrap_statistics(obs_real, rng, config.sim_spec.N_resamples);

            observable_mean_vector[k] = observable_mean;
            observable_std_vector[k] = observable_std;
            binder_mean_vector[k] = binder_mean;
            binder_std_vector[k] = binder_std;
            observable_autocorrelation_time_vector[k] 
            = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
        } else if (obs_type_vec[k] == "fredenhagen_marcu") {
            const auto& series = observable_vector[k];
            const size_t N = series.size();
            std::vector<double> obs_real, obs_imag;
            obs_real.reserve(N);
            obs_imag.reserve(N);

            for (auto const& v : series) {
                if (auto p = std::get_if<std::complex<double>>(&v)) {
                    obs_real.push_back(p->real());
                    obs_imag.push_back(p->imag());
                }
                else {
                    double d = std::get<double>(v);
                    obs_real.push_back(d);
                    obs_imag.push_back(0.0);
                }
            }

            const auto& [observable_mean, observable_std, binder_mean, binder_std] 
            = paratoric::statistics::get_bootstrap_statistics_fm(obs_real, obs_imag, rng, config.sim_spec.N_resamples);

            observable_mean_vector[k] = observable_mean;
            observable_std_vector[k] = observable_std;
            binder_mean_vector[k] = binder_mean;
            binder_std_vector[k] = binder_std;
            observable_autocorrelation_time_vector[k] 
            = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
        } else if (obs_type_vec[k] == "susceptibility") {
            const auto& series = observable_vector[k];
            const size_t N = series.size();
            std::vector<double> obs_real, obs_imag;
            obs_real.reserve(N);
            obs_imag.reserve(N);

            for (auto const& v : series) {
                if (auto p = std::get_if<std::complex<double>>(&v)) {
                    obs_real.push_back(p->real());
                    obs_imag.push_back(p->imag());
                }
                else {
                    double d = std::get<double>(v);
                    obs_real.push_back(d);
                    obs_imag.push_back(0.0);
                }
            }
            // TODO: Store susceptibility reducers in the observable registry.
            if ((config.sim_spec.observables[k] == "sigma_z_static_susceptibility" && Basis == 'x')) {
                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::bootstrap_offdiag_susceptibility(
                    obs_real, config.lat_spec.beta, config.param_spec.lmbda, 
                    lat.get_edge_count(), rng, config.sim_spec.N_resamples
                );
                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
            } else if (config.sim_spec.observables[k] == "sigma_x_static_susceptibility" && Basis == 'z') {
                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::bootstrap_offdiag_susceptibility(
                    obs_real, config.lat_spec.beta, config.param_spec.h, 
                    lat.get_edge_count(), rng, config.sim_spec.N_resamples
                );
                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
            } else if ((config.sim_spec.observables[k] == "sigma_z_dynamical_susceptibility" && Basis == 'x')
                        || (config.sim_spec.observables[k] == "sigma_x_dynamical_susceptibility" && Basis == 'z')) {
                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::bootstrap_offdiag_dynamical_susceptibility(
                    obs_real, obs_imag, lat.get_edge_count(), rng, config.sim_spec.N_resamples
                );
                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
            } else {
                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::get_bootstrap_statistics_susceptibility(
                    obs_real, obs_imag, rng, config.sim_spec.N_resamples
                );
                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(obs_real, config.sim_spec.observables[k]);
            }
        }
    }

    integrated_pot_energy_check = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );

    if (!almost_equal(integrated_pot_energy, integrated_pot_energy_check, 1e-5, 1e-13)) {
        throw std::runtime_error(std::format("Integrated potential energy mismatch. {} does not match {}.", integrated_pot_energy, integrated_pot_energy_check));
    }

    return Result{
        .series=std::move(observable_vector), 
        .mean=std::move(observable_mean_vector), 
        .mean_std=std::move(observable_std_vector), 
        .binder=std::move(binder_mean_vector), 
        .binder_std=std::move(binder_std_vector), 
        .tau_int=std::move(observable_autocorrelation_time_vector)
    };                                       
}

template<char Basis>
requires ValidBasis<Basis>
Result ExtendedToricCodeQMC<Basis>::get_hysteresis(
    const Config& config
) {

    detail::validate_hysteresis_config(config);

    if constexpr (Basis != 'x' && Basis != 'z') {
        throw std::invalid_argument("Basis must be either \"x\" or \"z\".");
    }

    if (Basis != config.lat_spec.basis) {
        throw std::invalid_argument("Template parameter basis and config.lat_spec basis must match.");
    }
    
    double h{}, lmbda{};
    
    if (!config.param_spec.h_hys.empty()) {
        h = config.param_spec.h_hys.front();
    } else {
        throw std::invalid_argument("h_hys must be non-empty.");
    }

    if (!config.param_spec.lmbda_hys.empty()) {
        lmbda = config.param_spec.lmbda_hys.front();
    } else {
        throw std::invalid_argument("lmbda_hys must be non-empty.");
    }

    if constexpr (Basis == 'x') {
        if (config.param_spec.J < 0) {
            throw std::invalid_argument("J must be non-negative in the x-basis.");
        } else if (lmbda < 0) {
            throw std::invalid_argument("lmbda must be non-negative in the x-basis.");
        }
    } else if constexpr (Basis == 'z') {
        if (config.param_spec.mu < 0) {
            throw std::invalid_argument("mu must be non-negative in the z-basis.");
        } else if (h < 0) {
            throw std::invalid_argument("h must be non-negative in the z-basis.");
        }
    }

    if (config.sim_spec.seed != 0) rng->set_seed(config.sim_spec.seed);

    auto obs_func_vec = get_obs_func_vec(config.sim_spec.observables);
    auto obs_type_vec = get_obs_type_vec(config.sim_spec.observables);
    
    // Hysteresis series: [schedule point][observable][sample].
    std::vector<std::vector<std::vector<std::variant< std::complex<double>, double>>>> hys_vector;
    std::vector<std::vector<double>> hys_mean,
                                     hys_mean_std,
                                     hys_binder,
                                     hys_binder_std,
                                     hys_autocorrelation_time;

    auto lat = Lattice(config.lat_spec, rng);
    
    double integrated_pot_energy = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );
    double acc_ratio = 1.;

    // Thermalization 
    for (int i = 0; i < config.sim_spec.N_thermalization; ++i) {
        metropolis_step(
            lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, 
            config.param_spec.h, config.param_spec.mu, config.param_spec.J, 
            config.param_spec.lmbda
        );
    }

    double integrated_pot_energy_check = total_integrated_pot_energy(
        lat, config.param_spec.h, config.param_spec.mu, config.param_spec.J, config.param_spec.lmbda
    );

    if (!almost_equal(integrated_pot_energy, integrated_pot_energy_check, 1e-5, 1e-13)) {
        throw std::runtime_error(std::format("Integrated potential energy mismatch. {} does not match {}.", integrated_pot_energy, integrated_pot_energy_check));
    }

    int metropolis_step_count = 0;
    int reset_potential_energy_count = static_cast<int>(lat.get_edge_count()*100000);

    for (size_t n = 0; n < config.param_spec.h_hys.size(); n++) {
        // Series: [observable][sample], in the requested observable order.
        const auto sample_count = static_cast<size_t>(std::max(config.sim_spec.N_samples, 0));
        std::vector<std::vector<std::variant< std::complex<double>, double>>> observable_vector;
        observable_vector.reserve(config.sim_spec.observables.size());
        std::vector<double> observable_mean_vector(config.sim_spec.observables.size(), 0.), 
                            observable_std_vector(config.sim_spec.observables.size(), 0.), 
                            binder_mean_vector(config.sim_spec.observables.size(), 0.), 
                            binder_std_vector(config.sim_spec.observables.size(), 0.), 
                            observable_autocorrelation_time_vector(config.sim_spec.observables.size(), 0.);
        for (const auto& obs_func : config.sim_spec.observables) {
            UNUSED(obs_func);
            observable_vector.emplace_back();
            observable_vector.back().reserve(sample_count);
        } 
        const auto path_out = n < config.out_spec.paths_out.size()
            ? config.out_spec.paths_out[n]
            : std::filesystem::path{};
        h = config.param_spec.h_hys[ n ];
        lmbda = config.param_spec.lmbda_hys[ n ];

        reinitialize_potential_energy(
            lat, integrated_pot_energy, h, config.param_spec.mu,
            config.param_spec.J, lmbda
        );

        if constexpr (Basis == 'x') {
            if (config.param_spec.J < 0) {
                throw std::invalid_argument("J must be non-negative in the x-basis.");
            } else if (lmbda < 0) {
                throw std::invalid_argument("lmbda must be non-negative in the x-basis.");
            }
        } else if constexpr (Basis == 'z') {
            if (config.param_spec.mu < 0) {
                throw std::invalid_argument("mu must be non-negative in the z-basis.");
            } else if (h < 0) {
                throw std::invalid_argument("h must be non-negative in the z-basis.");
            }
        }
        
        int N_rethermalization = static_cast<int>(config.sim_spec.N_thermalization/4);

        for (int t = 0; t < N_rethermalization; ++t) {
            ++metropolis_step_count;
            metropolis_step(
                lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, h, 
                config.param_spec.mu, config.param_spec.J, lmbda
            );
            if (metropolis_step_count == reset_potential_energy_count) [[unlikely]] {
                metropolis_step_count = 0;
                // Rebuild caches periodically to limit accumulated rounding error.
                reinitialize_potential_energy(
                    lat, integrated_pot_energy, h, config.param_spec.mu,
                    config.param_spec.J, lmbda
                );
            }
        }

        for (int i = 0; i < config.sim_spec.N_samples; ++i) {
            for (int j = 0; j < config.sim_spec.N_between_samples; ++j) {
                ++metropolis_step_count;
                metropolis_step(
                    lat, integrated_pot_energy, acc_ratio, config.lat_spec.beta, 
                    h, config.param_spec.mu, config.param_spec.J, lmbda
                );
                if (metropolis_step_count == reset_potential_energy_count) [[unlikely]] {
                    metropolis_step_count = 0;
                    // Rebuild caches periodically to limit accumulated rounding error.
                    reinitialize_potential_energy(
                        lat, integrated_pot_energy, h, config.param_spec.mu,
                        config.param_spec.J, lmbda
                    );
                }
            }

            for (size_t k = 0; k < config.sim_spec.observables.size(); k++)
                observable_vector[k].emplace_back(obs_func_vec[k](lat, h, lmbda, config.param_spec.mu, config.param_spec.J));

            if (config.out_spec.save_snapshots)
                lat.update_spin_string();
        }

        if (config.out_spec.save_snapshots) {
            lat.write_graph("snapshots", path_out);
        }

        for (size_t k = 0; k < config.sim_spec.observables.size(); k++) {
            if (obs_type_vec[k] == "real") {
                std::vector<double> obs_real;
                const auto& series = observable_vector[k];
                obs_real.reserve(series.size());
                for (auto const& x : series) {
                obs_real.emplace_back(std::get<double>(x));
                }

                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::get_bootstrap_statistics(obs_real, rng, config.sim_spec.N_resamples);

                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(
                    obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                );
            } else if (obs_type_vec[k] == "fredenhagen_marcu") {
                const auto& series = observable_vector[k];
                const size_t N = series.size();
                std::vector<double> obs_real, obs_imag;
                obs_real.reserve(N);
                obs_imag.reserve(N);

                for (auto const& v : series) {
                    if (auto p = std::get_if<std::complex<double>>(&v)) {
                        obs_real.push_back(p->real());
                        obs_imag.push_back(p->imag());
                    }
                    else {
                        double d = std::get<double>(v);
                        obs_real.push_back(d);
                        obs_imag.push_back(0.0);
                    }
                }

                const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                = paratoric::statistics::get_bootstrap_statistics_fm(obs_real, obs_imag, rng, config.sim_spec.N_resamples);

                observable_mean_vector[k] = observable_mean;
                observable_std_vector[k] = observable_std;
                binder_mean_vector[k] = binder_mean;
                binder_std_vector[k] = binder_std;
                observable_autocorrelation_time_vector[k] 
                = calculate_autocorrelation_time_with_warning(
                    obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                );
            } else if (obs_type_vec[k] == "susceptibility") {
                const auto& series = observable_vector[k];
                const size_t N = series.size();
                std::vector<double> obs_real, obs_imag;
                obs_real.reserve(N);
                obs_imag.reserve(N);

                for (auto const& v : series) {
                    if (auto p = std::get_if<std::complex<double>>(&v)) {
                        obs_real.push_back(p->real());
                        obs_imag.push_back(p->imag());
                    }
                    else {
                        double d = std::get<double>(v);
                        obs_real.push_back(d);
                        obs_imag.push_back(0.0);
                    }
                }

                if ((config.sim_spec.observables[k] == "sigma_z_static_susceptibility" && Basis == 'x')) {
                    const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                    = paratoric::statistics::bootstrap_offdiag_susceptibility(
                        obs_real, config.lat_spec.beta, lmbda, 
                        lat.get_edge_count(), rng, config.sim_spec.N_resamples
                    );
                    observable_mean_vector[k] = observable_mean;
                    observable_std_vector[k] = observable_std;
                    binder_mean_vector[k] = binder_mean;
                    binder_std_vector[k] = binder_std;
                    observable_autocorrelation_time_vector[k] 
                    = calculate_autocorrelation_time_with_warning(
                        obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                    );
                } else if (config.sim_spec.observables[k] == "sigma_x_static_susceptibility" && Basis == 'z') {
                    const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                    = paratoric::statistics::bootstrap_offdiag_susceptibility(
                        obs_real, config.lat_spec.beta, h, 
                        lat.get_edge_count(), rng, config.sim_spec.N_resamples
                    );
                    observable_mean_vector[k] = observable_mean;
                    observable_std_vector[k] = observable_std;
                    binder_mean_vector[k] = binder_mean;
                    binder_std_vector[k] = binder_std;
                    observable_autocorrelation_time_vector[k] 
                    = calculate_autocorrelation_time_with_warning(
                        obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                    );
                } else if ((config.sim_spec.observables[k] == "sigma_z_dynamical_susceptibility" && Basis == 'x')
                            || (config.sim_spec.observables[k] == "sigma_x_dynamical_susceptibility" && Basis == 'z')) {
                    const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                    = paratoric::statistics::bootstrap_offdiag_dynamical_susceptibility(
                        obs_real, obs_imag, lat.get_edge_count(), rng, config.sim_spec.N_resamples
                    );
                    observable_mean_vector[k] = observable_mean;
                    observable_std_vector[k] = observable_std;
                    binder_mean_vector[k] = binder_mean;
                    binder_std_vector[k] = binder_std;
                    observable_autocorrelation_time_vector[k] 
                    = calculate_autocorrelation_time_with_warning(
                        obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                    );
                } else {
                    const auto& [observable_mean, observable_std, binder_mean, binder_std] 
                    = paratoric::statistics::get_bootstrap_statistics_susceptibility(
                        obs_real, obs_imag, rng, config.sim_spec.N_resamples
                    );
                    observable_mean_vector[k] = observable_mean;
                    observable_std_vector[k] = observable_std;
                    binder_mean_vector[k] = binder_mean;
                    binder_std_vector[k] = binder_std;
                    observable_autocorrelation_time_vector[k] 
                    = calculate_autocorrelation_time_with_warning(
                        obs_real, config.sim_spec.observables[k], true, n, h, lmbda
                    );
                }
            } 
        }

        hys_vector.emplace_back( std::move(observable_vector) );
        hys_mean.emplace_back(std::move(observable_mean_vector));
        hys_mean_std.emplace_back(std::move(observable_std_vector));
        hys_binder.emplace_back(std::move(binder_mean_vector));
        hys_binder_std.emplace_back(std::move(binder_std_vector));
        hys_autocorrelation_time.emplace_back(std::move(observable_autocorrelation_time_vector));

        integrated_pot_energy_check = total_integrated_pot_energy(
            lat, h, config.param_spec.mu, config.param_spec.J, lmbda
        );

        if (!almost_equal(integrated_pot_energy, integrated_pot_energy_check, 1e-5, 1e-13)) {
            throw std::runtime_error(std::format("Integrated potential energy mismatch. {} does not match {}.", integrated_pot_energy, integrated_pot_energy_check));
        }
    }

    return Result{
        .series_hys=std::move(hys_vector), 
        .mean_hys=std::move(hys_mean), 
        .mean_std_hys=std::move(hys_mean_std), 
        .binder_hys=std::move(hys_binder), 
        .binder_std_hys=std::move(hys_binder_std), 
        .tau_int_hys=std::move(hys_autocorrelation_time)
    };
}

} // namespace paratoric
