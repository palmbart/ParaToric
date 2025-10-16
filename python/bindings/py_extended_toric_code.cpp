// ParaToric - Continuous-time QMC for the extended toric code in the x/z-basis
// Copyright (C) 2022-2025  Simon Mathias Linsel, Lode Pollet

#include "paratoric/mcmc/extended_toric_code.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include <complex>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <variant>
#include <optional>
#include <vector>

namespace py = pybind11;
namespace fs = std::filesystem;

using paratoric::ExtendedToricCode;

namespace {

inline std::complex<double>
to_complex(const std::variant<std::complex<double>, double>& v)
{
    if (auto p = std::get_if<std::complex<double>>(&v)) return *p;
    return std::complex<double>(std::get<double>(v), 0.0);
}

// 2D: [n_obs][n_samples] -> np.ndarray (n_obs, n_samples) complex128
py::array to_numpy_2d_complex(
    const std::vector<std::vector<std::variant<std::complex<double>, double>>>& vv)
{
    const ssize_t n_obs = static_cast<ssize_t>(vv.size());
    const ssize_t n_samples = n_obs ? static_cast<ssize_t>(vv[0].size()) : 0;

    for (const auto& row : vv)
        if (static_cast<ssize_t>(row.size()) != n_samples)
            throw std::runtime_error("observable series are not of equal length");

    py::array_t<std::complex<double>> a({n_obs, n_samples});
    auto buf = a.mutable_unchecked<2>();  // typed array -> OK in pybind11 v3
    for (ssize_t i = 0; i < n_obs; ++i)
        for (ssize_t j = 0; j < n_samples; ++j)
            buf(i, j) = to_complex(vv[i][static_cast<size_t>(j)]);
    return a;
}

// 3D: [n_steps][n_obs][n_samples] -> np.ndarray (n_steps, n_obs, n_samples) complex128
py::array to_numpy_3d_complex(
    const std::vector<std::vector<std::vector<std::variant<std::complex<double>, double>>>>& vvv)
{
    const ssize_t n_steps   = static_cast<ssize_t>(vvv.size());
    const ssize_t n_obs     = n_steps ? static_cast<ssize_t>(vvv[0].size()) : 0;
    const ssize_t n_samples = (n_steps && n_obs) ? static_cast<ssize_t>(vvv[0][0].size()) : 0;

    for (const auto& mat : vvv) {
        if (static_cast<ssize_t>(mat.size()) != n_obs)
            throw std::runtime_error("hysteresis: each step must have same number of observables");
        for (const auto& row : mat)
            if (static_cast<ssize_t>(row.size()) != n_samples)
                throw std::runtime_error("hysteresis: each observable series must have same length");
    }

    py::array_t<std::complex<double>> a({n_steps, n_obs, n_samples});
    auto buf = a.mutable_unchecked<3>();
    for (ssize_t s = 0; s < n_steps; ++s)
        for (ssize_t i = 0; i < n_obs; ++i)
            for (ssize_t j = 0; j < n_samples; ++j)
                buf(s, i, j) = to_complex(vvv[static_cast<size_t>(s)][static_cast<size_t>(i)][static_cast<size_t>(j)]);
    return a;
}

// 1D vector<double> -> np.ndarray (n,)
py::array to_numpy_1d_double(const std::vector<double>& v)
{
    const ssize_t n = static_cast<ssize_t>(v.size());
    py::array_t<double> a(py::array::ShapeContainer{static_cast<py::ssize_t>(n)});
    auto buf = a.mutable_unchecked<1>();
    for (ssize_t i = 0; i < n; ++i) buf(i) = v[static_cast<size_t>(i)];
    return a;
}


// 2D vector<vector<double>> -> np.ndarray (n_rows, n_cols)
py::array to_numpy_2d_double(const std::vector<std::vector<double>>& vv)
{
    const ssize_t n_rows = static_cast<ssize_t>(vv.size());
    const ssize_t n_cols = n_rows ? static_cast<ssize_t>(vv[0].size()) : 0;

    for (const auto& row : vv)
        if (static_cast<ssize_t>(row.size()) != n_cols)
            throw std::runtime_error("matrix rows have different lengths");

    py::array_t<double> a({n_rows, n_cols});
    auto buf = a.mutable_unchecked<2>();
    for (ssize_t i = 0; i < n_rows; ++i)
        for (ssize_t j = 0; j < n_cols; ++j)
            buf(i, j) = vv[static_cast<size_t>(i)][static_cast<size_t>(j)];
    return a;
}

} // namespace

PYBIND11_MODULE(_paratoric, m) {
    m.doc() = R"pbdoc(
ParaToric bindings that return NumPy arrays.

Functions
---------
get_thermalization
get_sample
get_hysteresis
)pbdoc";

    auto etc = m.def_submodule(
        "extended_toric_code",
        "Extended toric code QMC"
    );

    // get_thermalization -> (series[n_obs, n_therm], acc_ratio[n_therm])
    etc.def(
        "get_thermalization",
        [](int N_thermalization,
           int N_resamples,
           const std::vector<std::string>& observables,
           int seed,
           double mu,
           double h,
           double J,
           double lmbda,
           char basis,
           const std::string& lattice_type,
           int system_size,
           double beta,
           const std::string& boundaries,
           int default_spin,
           bool save_snapshots,
           std::optional<fs::path> path_out_opt) {
            if (basis != 'x' && basis != 'z')
                throw std::invalid_argument("Basis must be 'x' or 'z'.");
            fs::path path_out = path_out_opt.value_or(fs::path{});

            std::vector<std::vector<std::variant<std::complex<double>, double>>> series;
            std::vector<double> acc_ratio;
            {
                py::gil_scoped_release release;
                paratoric::SimSpec sim_spec = {
                    .N_thermalization = N_thermalization,
                    .N_resamples = N_resamples,
                    .seed = seed,
                    .observables = observables
                };
                paratoric::ParamSpec param_spec = {
                    .mu = mu,
                    .h = h,
                    .J = J,
                    .lmbda = lmbda
                };
                paratoric::LatSpec lat_spec = {
                    basis,
                    lattice_type,
                    system_size,
                    beta,
                    boundaries,
                    default_spin
                };
                paratoric::OutSpec out_spec = {
                    .path_out = path_out,
                    .save_snapshots = save_snapshots
                };

                auto result = ExtendedToricCode::get_thermalization(
                    paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
                series = std::move(result.series);
                acc_ratio = std::move(result.acc_ratio);
            }
            py::tuple out(2);
            out[0] = to_numpy_2d_complex(series);
            out[1] = to_numpy_1d_double(acc_ratio);
            return out;
        },
        py::arg("N_thermalization"),
        py::arg("N_resamples") = 1000,
        py::arg("observables"),
        py::arg("seed") = 0,
        py::arg("mu"),
        py::arg("h"),
        py::arg("J"),
        py::arg("lmbda"),
        py::arg("basis") = 'x',
        py::arg("lattice_type"),
        py::arg("system_size"),
        py::arg("beta"),
        py::arg("boundaries") = "periodic",
        py::arg("default_spin") = 1,
        py::arg("save_snapshots") = false,
        py::arg("path_out") = py::none(),
        py::doc(R"pbdoc(
Run a QMC thermalization of the extended toric code and return observables and acceptance ratio diagnostics.

Parameters
----------
N_thermalization : int
    Number of thermalization steps.
N_resamples : int, optional
    Number of bootstrap resamples.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
mu : float
    Hamiltonian parameter (star term).
h : float
    Hamiltonian parameter (electric field term).
J : float
    Hamiltonian parameter (plaquette term).
lmbda : float
    Hamiltonian parameter (gauge field term).
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
beta : float
    Inverse temperature.
boundaries : str, optional
    Boundary condition ("periodic", "open").
default_spin : int, optional
    Default link spin (+1 or -1).
save_snapshots : bool, optional
    If True, snapshots may be saved (not required for array returns).
path_out : pathlib.Path | None, optional
    Output directory for snapshots when saving; ignored for in-memory arrays.

Returns
-------
series : numpy.ndarray (complex128), shape (n_obs, N_thermalization)
    Time series per observable during thermalization.
acc_ratio : numpy.ndarray (float64), shape (N_thermalization,)
    Acceptance ratio diagnostics during thermalization.

Notes
-----
)pbdoc")
    );

    // get_sample -> (series[n_obs, n_samples], mean[n_obs], std[n_obs], binder[n_obs], binder_std[n_obs], tau_int[n_obs])
    etc.def(
        "get_sample",
        [](int N_samples,
           int N_thermalization,
           int N_between_samples,
           int N_resamples,
           bool custom_therm,
           const std::vector<std::string>& observables,
           int seed,
           double mu,
           double h,
           double h_therm,
           double J,
           double lmbda,
           double lmbda_therm,
           char basis,
           const std::string& lattice_type,
           int system_size,
           double beta,
           const std::string& boundaries,
           int default_spin,
           bool save_snapshots,
           std::optional<fs::path> path_out_opt) {
            if (basis != 'x' && basis != 'z')
                throw std::invalid_argument("basis must be 'x' or 'z'");
            fs::path path_out = path_out_opt.value_or(fs::path{});

            std::vector<std::vector<std::variant<std::complex<double>, double>>> series;
            std::vector<double> mean, mean_std, binder, binder_std, tau_int;
            {
                py::gil_scoped_release release;
                paratoric::SimSpec sim_spec = {
                    .N_samples = N_samples,
                    .N_thermalization = N_thermalization,
                    .N_between_samples = N_between_samples,
                    .N_resamples = N_resamples,
                    .custom_therm = custom_therm,
                    .seed = seed,
                    .observables = observables
                };
                paratoric::ParamSpec param_spec = {
                    .mu = mu,
                    .h = h,
                    .J = J,
                    .lmbda = lmbda,
                    .h_therm = h_therm,
                    .lmbda_therm = lmbda_therm
                };
                paratoric::LatSpec lat_spec = {
                    basis,
                    lattice_type,
                    system_size,
                    beta,
                    boundaries,
                    default_spin
                };
                paratoric::OutSpec out_spec = {
                    .path_out = path_out,
                    .save_snapshots = save_snapshots
                };

                auto result = ExtendedToricCode::get_sample(
                    paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
                series      = std::move(result.series);
                mean        = std::move(result.mean);
                mean_std    = std::move(result.mean_std);
                binder      = std::move(result.binder);
                binder_std  = std::move(result.binder_std);
                tau_int     = std::move(result.tau_int);
            }

            py::tuple out(6);
            out[0] = to_numpy_2d_complex(series);
            out[1] = to_numpy_1d_double(mean);
            out[2] = to_numpy_1d_double(mean_std);
            out[3] = to_numpy_1d_double(binder);
            out[4] = to_numpy_1d_double(binder_std);
            out[5] = to_numpy_1d_double(tau_int);
            return out;
        },
        py::arg("N_samples"),
        py::arg("N_thermalization"),
        py::arg("N_between_samples"),
        py::arg("N_resamples") = 1000,
        py::arg("custom_therm") = false,
        py::arg("observables"),
        py::arg("seed") = 0,
        py::arg("mu"),
        py::arg("h"),
        py::arg("h_therm") = 0,
        py::arg("J"),
        py::arg("lmbda"),
        py::arg("lmbda_therm") = 0,
        py::arg("basis") = 'x',
        py::arg("lattice_type"),
        py::arg("system_size"),
        py::arg("beta"),
        py::arg("boundaries") = "periodic",
        py::arg("default_spin") = 1,
        py::arg("save_snapshots") = false,
        py::arg("path_out") = py::none(),
        py::doc(R"pbdoc(
Run a QMC simulation of the extended toric code and return observables and statistics.

Parameters
----------
N_samples : int
    Number of snapshots (stored samples).
N_thermalization : int
    Number of thermalization steps before sampling.
N_between_samples : int
    Steps between stored samples (decorrelation).
N_resamples : int, optional
    Number of bootstrap resamples.
custom_therm : bool, optional
    If True, use (h_therm, lmbda_therm) during thermalization.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
mu : float
    Hamiltonian parameter (star term).
h : float
    Hamiltonian parameter (electric field term).
h_therm : float
    Electric field used during thermalization if `custom_therm=True`.
J : float
    Hamiltonian parameter (plaquette term).
lmbda : float
    Hamiltonian parameter (gauge field term).
lmbda_therm : float
    Gauge field used during thermalization if `custom_therm=True`.
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
beta : float
    Inverse temperature.
boundaries : str, optional
    Boundary condition ("periodic", "open").
default_spin : int, optional
    Default link spin (+1 or -1).
save_snapshots : bool, optional
    If True, snapshots may be saved (not required for array returns).
path_out : pathlib.Path | None, optional
    Output directory for snapshots when saving; ignored for in-memory arrays.

Returns
-------
series : numpy.ndarray (complex128), shape (n_obs, N_samples)
    Time series per observable.
mean : numpy.ndarray (float64), shape (n_obs,)
    Sample mean per observable.
std : numpy.ndarray (float64), shape (n_obs,)
    Bootstrap standard deviation per observable.
binder : numpy.ndarray (float64), shape (n_obs,)
    Binder ratio (or analogous) per observable.
binder_std : numpy.ndarray (float64), shape (n_obs,)
    Bootstrap error of the Binder ratio.
tau_int : numpy.ndarray (float64), shape (n_obs,)
    Integrated autocorrelation time estimate.

Notes
-----
Complex observables are returned as complex128; real ones have zero imaginary part.
)pbdoc")
    );

    // get_hysteresis -> (
    //   series3d[n_steps, n_obs, n_samples],
    //   mean2d[n_steps, n_obs], std2d[n_steps, n_obs],
    //   binder2d[n_steps, n_obs], binder_std2d[n_steps, n_obs],
    //   tau2d[n_steps, n_obs]
    // )
    etc.def(
        "get_hysteresis",
        [](int N_samples,
           int N_thermalization,
           int N_between_samples,
           int N_resamples,
           const std::vector<std::string>& observables,
           int seed,
           double mu,
           const std::vector<double>& h_hys,
           double J,
           const std::vector<double>& lmbda_hys,
           char basis,
           const std::string& lattice_type,
           int system_size,
           double beta,
           const std::string& boundaries,
           int default_spin,
           bool save_snapshots,
           std::optional<std::vector<fs::path>> paths_out_opt) {
            if (basis != 'x' && basis != 'z')
                throw std::invalid_argument("basis must be 'x' or 'z'");

            std::vector<fs::path> paths_out;
            if (paths_out_opt.has_value()) {
                paths_out = std::move(*paths_out_opt);
            } else {
                paths_out.resize(h_hys.size()); // empty paths
            }
            if (save_snapshots && paths_out.size() != h_hys.size())
                throw std::invalid_argument("paths_out must match h_hys size when save_snapshots=True");

            std::vector<std::vector<std::vector<std::variant<std::complex<double>, double>>>> series3d;
            std::vector<std::vector<double>> mean2d, std2d, binder2d, binder_std2d, tau2d;

            {
                py::gil_scoped_release release;
                paratoric::SimSpec sim_spec = {
                    .N_samples = N_samples,
                    .N_thermalization = N_thermalization,
                    .N_between_samples = N_between_samples,
                    .N_resamples = N_resamples,
                    .seed = seed,
                    .observables = observables
                };
                paratoric::ParamSpec param_spec = {
                    .mu = mu,
                    .J = J,
                    .h_hys = h_hys,
                    .lmbda_hys = lmbda_hys
                };
                paratoric::LatSpec lat_spec = {
                    basis,
                    lattice_type,
                    system_size,
                    beta,
                    boundaries,
                    default_spin
                };
                paratoric::OutSpec out_spec = {
                    .paths_out = paths_out,
                    .save_snapshots = save_snapshots
                };

                auto result = ExtendedToricCode::get_hysteresis(
                    paratoric::Config{sim_spec, param_spec, lat_spec, out_spec}
                );
                series3d        = std::move(result.series_hys);
                mean2d          = std::move(result.mean_hys);
                std2d           = std::move(result.mean_std_hys);
                binder2d        = std::move(result.binder_hys);
                binder_std2d    = std::move(result.binder_std_hys);
                tau2d           = std::move(result.tau_int_hys);
            }

            py::tuple out(6);
            out[0] = to_numpy_3d_complex(series3d);
            out[1] = to_numpy_2d_double(mean2d);
            out[2] = to_numpy_2d_double(std2d);
            out[3] = to_numpy_2d_double(binder2d);
            out[4] = to_numpy_2d_double(binder_std2d);
            out[5] = to_numpy_2d_double(tau2d);
            return out;
        },
        py::arg("N_samples"),
        py::arg("N_thermalization"),
        py::arg("N_between_samples"),
        py::arg("N_resamples") = 1000,
        py::arg("observables"),
        py::arg("seed") = 0,
        py::arg("mu"),
        py::arg("h_hys"),
        py::arg("J"),
        py::arg("lmbda_hys"),
        py::arg("basis") = 'x',
        py::arg("lattice_type"),
        py::arg("system_size"),
        py::arg("beta"),
        py::arg("boundaries") = "periodic",
        py::arg("default_spin") = 1,
        py::arg("save_snapshots") = false,
        py::arg("paths_out") = py::none(),
        py::doc(R"pbdoc(
Run a QMC hysteresis sweep of the extended toric code and return stacked observables across steps.

Parameters
----------
N_samples : int
    Number of snapshots per hysteresis step.
N_thermalization : int
    Number of thermalization steps before sampling.
N_between_samples : int
    Steps between stored samples (decorrelation).
N_resamples : int, optional
    Number of bootstrap resamples.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
mu : float
    Hamiltonian parameter (star term).
h_hys : list[float]
    Electric field values for each hysteresis step.
J : float
    Hamiltonian parameter (plaquette term).
lmbda_hys : list[float]
    Gauge field values for each hysteresis step.
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
beta : float
    Inverse temperature.
boundaries : str, optional
    Boundary condition ("periodic", "open").
default_spin : int, optional
    Default link spin (+1 or -1).
save_snapshots : bool, optional
    If True, snapshots may be saved (not required for array returns).
paths_out : list[pathlib.Path] | None, optional
    Output directories per step when saving; ignored for in-memory arrays.

Returns
-------
series3d : numpy.ndarray (complex128), shape (n_steps, n_obs, N_samples)
    Time series per observable for each hysteresis step.
mean2d : numpy.ndarray (float64), shape (n_steps, n_obs)
std2d : numpy.ndarray (float64), shape (n_steps, n_obs)
binder2d : numpy.ndarray (float64), shape (n_steps, n_obs)
binder_std2d : numpy.ndarray (float64), shape (n_steps, n_obs)
tau2d : numpy.ndarray (float64), shape (n_steps, n_obs)

Notes
-----
The number of steps `n_steps` equals `len(h_hys)` (and `len(lmbda_hys)`).
)pbdoc")
    );
}
