// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    PTC_STATUS_OK = 0,
    PTC_STATUS_INVALID_ARGUMENT = 1,
    PTC_STATUS_RUNTIME_ERROR = 2,
    PTC_STATUS_NO_MEMORY = 3,
    PTC_STATUS_INTERNAL_ERROR = 4
} ptc_status_t;

/**
 * @brief Get a thread-local, nul-terminated error string for the last error.
 * 
 * @return error string
 */
const char* ptc_last_error(void);


/**
 * @brief complex number
 * 
 * @param re real part
 * @param im imaginary part
 * 
 */
typedef struct { double re, im; } ptc_complex_t;

/**
 * @brief Either real (0) or complex (1)
 * 
 */
typedef enum {
    PTC_SCALAR_REAL   = 0,
    PTC_SCALAR_COMPLEX = 1
} ptc_scalar_kind_t;

/**
 * @brief Can either hold real or complex number 
 * 
 * @param kind either real (0) or complex (1)
 * @param re real part
 * @param im imaginary part
 * 
 */
typedef struct {
    ptc_scalar_kind_t kind;
    /* When kind==REAL, use 're' only; when COMPLEX use 're' and 'im'. */
    double re, im;
} ptc_scalar_t;

/**
 * @brief double vector
 * 
 * @param data holds double values
 * @param len number of elements in data
 */
typedef struct {
    double* data;
    size_t  len;
} ptc_dvec_t;

/**
 * @brief ptc_scalar_t vector
 * 
 * @param data holds ptc_scalar_t values
 * @param len number of elements in data
 * 
 */
typedef struct {
    ptc_scalar_t* data;
    size_t        len;
} ptc_scalar_vec_t;

/**
 * @brief Vector of ptc_scalar_t vectors
 * 
 * @param rows holds ptc_scalar_t vectors
 * @param nrows number of elements in rows
 * 
 */
typedef struct {
    ptc_scalar_vec_t* rows; /* rows[i] is one time series / operator series */
    size_t            nrows;
} ptc_series_t;

/**
 * @brief Vector of ptc_series_t vectors
 * 
 * @param blocks holds ptc_series_t vectors
 * @param nblocks number of elements in blocks
 * 
 */
typedef struct {
    ptc_series_t* blocks;
    size_t        nblocks;
} ptc_series_blocks_t;

/**
 * @brief Vector of ptc_dvec_t vectors
 * 
 * @param rows holds ptc_dvec_t vectors
 * @param nrows number of elements in rows
 * 
 */
typedef struct {
    ptc_dvec_t* rows;
    size_t      nrows;
} ptc_dmat_t;

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
typedef struct {
    char  basis;            
    const char* lattice_type; 
    int   system_size;
    double beta;
    const char* boundaries; 
    int   default_spin;
} ptc_lat_spec_t;

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
 * @param h_hys_len number of elements in h_hys
 * @param lmbda_hys the Hamiltonian hysteresis parameters (gauge field term)
 * @param lmbda_hys_len number of elements in lmbda_hys
 * 
 */
typedef struct {
    double mu;
    double h;
    double J;
    double lmbda;

    /* Optional “thermalization schedule” values. Use NaN if unset. */
    double h_therm;
    double lmbda_therm;

    /* Optional hysteresis schedules. Null or len=0 if unused. */
    const double* h_hys;     size_t h_hys_len;
    const double* lmbda_hys; size_t lmbda_hys_len;
} ptc_param_spec_t;

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
 * @param N_observables number of elements in observables
 * 
 */
typedef struct {
    int  N_samples;         /* 0 = thermalization-only */
    int  N_thermalization;
    int  N_between_samples;
    int  N_resamples;
    bool custom_therm;      /* used only by sampling */
    int  seed;              /* 0 = random seed */

    /* Array of observable names; can be NULL/0. */
    const char* const* observables;
    size_t             N_observables;
} ptc_sim_spec_t;

/**
 * @brief Output specification
 * 
 * @param path_out the directory where output files will be written to 
 * @param paths_out the directories (hysteresis simulation) where output files will be written to
 * @param N_paths_out number of paths in paths_out
 * @param save_snapshots whether snapshots should be saved
 * 
 */
typedef struct {
    /* Single path or NULL. */
    const char* path_out;

    /* Multiple paths (for hysteresis) or NULL/0. */
    const char* const* paths_out;
    size_t             N_paths_out;

    bool save_snapshots;
} ptc_out_spec_t;

/**
 * @brief Configuration for starting a simulation
 * 
 * @param ptc_sim_spec_t Simulation parameters (N_samples, seed, ...) 
 * @param ptc_param_spec_t Hamiltonian parameters
 * @param ptc_lat_spec_t Lattice parameters (system_size, beta, ...)
 * @param ptc_out_spec_t Output parameters (path_out, ...)
 * 
 */
typedef struct {
    ptc_sim_spec_t     sim;
    ptc_param_spec_t   params;
    ptc_lat_spec_t     lat;
    ptc_out_spec_t     out;
} ptc_config_t;

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
typedef struct {
    /* Per-observable (outer), per-time-step (inner) series with scalar entries. */
    ptc_series_t series;

    /* Acceptance ratios per sample (if applicable); else can be empty. */
    ptc_dvec_t acc_ratio;

    /* Scalar summaries per observable (same length across these vectors). */
    ptc_dvec_t mean, mean_std, binder, binder_std, tau_int;

    /* Hysteresis versions (outer index = hysteresis step). */
    ptc_series_blocks_t series_hys;
    ptc_dmat_t mean_hys, mean_std_hys, binder_hys, binder_std_hys, tau_int_hys;
} ptc_result_t;

/**
 * @brief Release heap memory inside a result and zero it. Safe to call on a zeroed struct.
 * 
 * @param r the result structure 
 */
void ptc_result_destroy(ptc_result_t* r);


typedef struct ptc_handle_t ptc_handle_t;

/**
 * @brief Create the facade
 * 
 * @param out_handle pointer to ptc_handle_t, can be null pointer 
 * @return ptc_status_t 
 */
ptc_status_t ptc_create(ptc_handle_t** out_handle);

/**
 * @brief Destroy the facade
 * 
 * @param handle pointer to ptc_handle_t
 */
void         ptc_destroy(ptc_handle_t* handle);

/* ========= Main entry points (mirror the C++ façade methods) =========
 * Each fills 'out' on success. Caller owns ptc_result_destroy(out).
 */

/**
 * @brief Run extended toric code thermalization.
 * 
 * @param handle handle struct
 * 
 * @param config contains all parameters of the simulation
 * @param config.sim.N_thermalization the number of thermalization steps
 * @param config.sim.N_resamples the number of bootstrap resamples
 * @param config.sim.observables all observables that are calculated for every snapshot
 * @param config.sim.observables.N_observables number of elements in observables
 * @param config.sim.seed the seed for the pseudorandom number generator
 * @param config.params.mu the Hamiltonian parameter (star term)
 * @param config.params.h the Hamiltonian parameter (electric field term)
 * @param config.params.J the Hamiltonian parameter (plaquette term)
 * @param config.params.lmbda the Hamiltonian parameter (gauge field term)
 * @param config.lat.basis eigenbasis of the spins, either 'x' or 'z'
 * @param config.lat.lattice_type the lattice type, e.g. "triangular"
 * @param config.lat.system_size the system size of the lattice (in one dimension)
 * @param config.lat.beta inverse temperature
 * @param config.lat.boundaries the boundary condition of the lattice (periodic, open)
 * @param config.lat.default_spin the default spin on the links (1 or -1)
 * @param config.out.path_out output directory for snapshots 
 * @param config.out.save_snapshots whether snapshots should be saved (every 10000th snapshot will be saved)
 * 
 * @param out Filled on success.
 * @param out.series      Time series (per snapshot) of all requested observables,
 *                        measured during thermalization. Each observable’s entry
 *                        contains one value per recorded snapshot, in time order.
 * @param out.acc_ratio   Time series of Monte Carlo acceptance ratios.
 * 
 * @return ptc_status_t exit status 
 */
ptc_status_t ptc_get_thermalization(
    ptc_handle_t* handle,
    const ptc_config_t* config,
    ptc_result_t* out
);

/**
 * @brief Run extended toric code simulation.
 * 
 * @param handle handle struct
 * 
 * @param config contains all parameters of the simulation
 * @param config the configuration object
 * @param config.sim.N_samples the number of snapshots
 * @param config.sim.N_thermalization the number of thermalization steps
 * @param config.sim.N_between_samples the number of steps between snapshots
 * @param config.sim.N_resamples the number of bootstrap resamples
 * @param config.sim.custom_therm if custom thermalization is used (to probe hysteresis)
 * @param config.sim.observables all observables that are calculated for every snapshot
 * @param config.sim.observables.N_observables number of elements in observables
 * @param config.sim.seed the seed for the pseudorandom number generator
 * @param config.params.mu the Hamiltonian parameter (star term)
 * @param config.params.h the Hamiltonian parameter (electric field term)
 * @param config.params.J the Hamiltonian parameter (plaquette term)
 * @param config.params.lmbda the Hamiltonian parameter (gauge field term)
 * @param config.params.h_therm the thermalization value of h, used if custom_therm enabled
 * @param config.params.lmbda_therm the thermalization value of lmbda, used if custom_therm enabled
 * @param config.lat.basis eigenbasis of the spins, either 'x' or 'z'
 * @param config.lat.lattice_type the lattice type, e.g. "triangular"
 * @param config.lat.system_size the system size of the lattice (in one dimension)
 * @param config.lat.beta inverse temperature
 * @param config.lat.boundaries the boundary condition of the lattice (periodic, open)
 * @param config.lat.default_spin the default spin on the links (1 or -1)
 * @param config.out.path_out output directory for snapshots 
 * @param config.out.save_snapshots whether snapshots should be saved (every snapshot will be saved)
 * 
 * @param out Filled on success.
 * @param out.series      Full counting statistics of all requested observables in the input order,
 *                        Each observable’s entry contains one value per recorded snapshot, 
 *                        in time order.
 * @param out.mean        Bootstrap observable means
 * @param out.mean_std    Bootstrap standard errors of the mean
 * @param out.binder      Bootstrap binder ratios
 * @param out.binder_std  Bootstrap standard errors of the binder ratios
 * @param out.tau_int     Estimated integrated autocorrelation times
 * 
 * @return ptc_status_t  
 */
ptc_status_t ptc_get_sample(
    ptc_handle_t* handle,
    const ptc_config_t* config,
    ptc_result_t* out
);

/**
 * @brief Run extended toric code hysteresis simulation.
 * 
 * @param handle handle struct
 * 
 * @param config contains all parameters of the simulation
 * @param config the configuration object
 * @param config.sim.N_samples the number of snapshots
 * @param config.sim.N_thermalization the number of thermalization steps
 * @param config.sim.N_between_samples the number of steps between snapshots
 * @param config.sim.N_resamples the number of bootstrap resamples
 * @param config.sim.observables all observables that are calculated for every snapshot
 * @param config.sim.observables.N_observables number of elements in observables
 * @param config.sim.seed the seed for the pseudorandom number generator
 * @param config.params.mu the Hamiltonian parameter (star term)
 * @param config.params.h_hys the Hamiltonian parameters (electric field term)
 * @param config.params.J the Hamiltonian parameter (plaquette term)
 * @param config.params.lmbda_hys the Hamiltonian parameters (gauge field term)
 * @param config.lat.basis eigenbasis of the spins, either 'x' or 'z'
 * @param config.lat.lattice_type the lattice type, e.g. "triangular"
 * @param config.lat.system_size the system size of the lattice (in one dimension)
 * @param config.lat.beta inverse temperature
 * @param config.lat.boundaries the boundary condition of the lattice (periodic, open)
 * @param config.lat.default_spin the default spin on the links (1 or -1)
 * @param config.out.paths_out all output directories for snapshots 
 * @param config.out.N_paths_out number of elements in paths_out
 * @param config.out.save_snapshots whether snapshots should be saved (every snapshot will be saved)
 * 
 * @param out Filled on success.
 * @return out.series_hys      Full counting statistics of all requested observables.
 *                             Each element in the outer vector represents one parameter point (h_hys, lmbda_hys)
 *                             Each element in the middle vectors represent one observables in the order of input.
 *                             The inner vector are the full counting statistics in time order.
 * @return out.mean_hys        Bootstrap observable means. 
 *                             Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                             Each element in the inner vectors represents one observables in the order of input.
 * @return out.mean_std_hys    Bootstrap standard errors of the mean.
 *                             Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                             Each element in the inner vectors represents one observables in the order of input.
 * @return out.binder_hys      Bootstrap binder ratios.
 *                             Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                             Each element in the inner vectors represents one observables in the order of input.
 * @return out.binder_std_hys  Bootstrap standard errors of the binder ratios. 
 *                             Each element of the outer vector represents one p rameter point (h_hys, lmbda_hys).
 *                             Each element in the inner vectors represents one observables in the order of input.
 * @return out.tau_int_hys     Estimated integrated autocorrelation times.
 *                             Each element of the outer vector represents one parameter point (h_hys, lmbda_hys).
 *                             Each element in the inner vectors represents one observables in the order of input.
 * 
 * @return ptc_status_t 
 */
ptc_status_t ptc_get_hysteresis(
    ptc_handle_t* handle,
    const ptc_config_t* config,
    ptc_result_t* out
);

#ifdef __cplusplus
} /* extern "C" */
#endif
