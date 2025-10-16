# python/paratoric/_paratoric.pyi
from __future__ import annotations
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Sequence, Tuple, Optional

class extended_toric_code:
    def get_thermalization(
        N_thermalization: int,
        beta: float,
        mu: float,
        h: float,
        J: float,
        lmbda: float,
        N_resamples: int,
        observables: Sequence[str],
        seed: int = ...,
        basis: str = ...,
        lattice_type: str = ...,
        system_size: int = ...,
        boundaries: str = ...,
        default_spin: int = ...,
        save_snapshots: bool = ...,
        path_out: Optional[Path] = ...,
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
        """
    Run a QMC thermalization of the extended toric code and return observables and acceptance ratio diagnostics.

Parameters
----------
N_thermalization : int
    Number of thermalization steps.
beta : float
    Inverse temperature.
mu : float
    Hamiltonian parameter (star term).
h : float
    Hamiltonian parameter (electric field term).
J : float
    Hamiltonian parameter (plaquette term).
lmbda : float
    Hamiltonian parameter (gauge field term).
N_resamples : int
    Number of bootstrap resamples.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
boundaries : str
    Boundary condition ("periodic" or "open").
default_spin : int
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
-----"""

    def get_sample(
        N_samples: int,
        N_thermalization: int,
        N_between_samples: int,
        beta: float,
        mu: float,
        h: float,
        h_therm: float,
        J: float,
        lmbda: float,
        lmbda_therm: float,
        N_resamples: int,
        custom_therm: bool,
        observables: Sequence[str],
        seed: int = ...,
        basis: str = ...,
        lattice_type: str = ...,
        system_size: int = ...,
        boundaries: str = ...,
        default_spin: int = ...,
        save_snapshots: bool = ...,
        path_out: Optional[Path] = ...,
    ) -> tuple[
        npt.NDArray[np.complex128],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Run a QMC simulation of the extended toric code and return observables and statistics.

Parameters
----------
N_samples : int
    Number of snapshots (stored samples).
N_thermalization : int
    Number of thermalization steps before sampling.
N_between_samples : int
    Steps between stored samples (decorrelation).
beta : float
    Inverse temperature.
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
N_resamples : int
    Number of bootstrap resamples.
custom_therm : bool, optional
    If True, use (h_therm, lmbda_therm) during thermalization.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
boundaries : str
    Boundary condition ("periodic" or "open").
default_spin : int
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
Complex observables are returned as complex128; real ones have zero imaginary part."""

    def get_hysteresis(
        N_samples: int,
        N_thermalization: int,
        N_between_samples: int,
        beta: float,
        mu: float,
        h_hys: Sequence[float],
        J: float,
        lmbda_hys: Sequence[float],
        N_resamples: int,
        observables: Sequence[str],
        seed: int = ...,
        basis: str = ...,
        lattice_type: str = ...,
        system_size: int = ...,
        boundaries: str = ...,
        default_spin: int = ...,
        save_snapshots: bool = ...,
        paths_out: Optional[Sequence[Path]] = ...,
    ) -> tuple[
        npt.NDArray[np.complex128],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Run a QMC hysteresis sweep of the extended toric code and return stacked observables across steps.

Parameters
----------
N_samples : int
    Number of snapshots per hysteresis step.
N_thermalization : int
    Number of thermalization steps before sampling.
N_between_samples : int
    Steps between stored samples (decorrelation).
beta : float
    Inverse temperature.
mu : float
    Hamiltonian parameter (star term).
h_hys : list[float]
    Electric field values for each hysteresis step.
J : float
    Hamiltonian parameter (plaquette term).
lmbda_hys : list[float]
    Gauge field values for each hysteresis step.
N_resamples : int
    Number of bootstrap resamples.
observables : list[str]
    Observable names computed at every snapshot.
seed : int, optional
    PRNG seed. 0 will lead to random seed.
basis : {'x','z'}, optional
    Basis for the templated backend; selects the implementation at runtime.
lattice_type : str
    Lattice type (e.g. "triangular").
system_size : int
    Linear system size of the lattice.
boundaries : str
    Boundary condition ("periodic" or "open").
default_spin : int
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
The number of steps `n_steps` equals `len(h_hys)` (and `len(lmbda_hys)`)."""