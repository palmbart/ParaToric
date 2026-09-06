"""Typed array-returning interface; docstrings match the compiled bindings."""
from __future__ import annotations
import numpy as np
import numpy.typing as npt
from pathlib import Path
from typing import Sequence, Tuple, Optional

class extended_toric_code:
    """Type-checking namespace for the extension submodule."""

    @staticmethod
    def get_thermalization(
        N_thermalization: int,
        *,
        N_resamples: int = ...,
        observables: Sequence[str],
        seed: int = ...,
        mu: float,
        h: float,
        J: float,
        lmbda: float,
        basis: str = ...,
        lattice_type: str,
        system_size: int,
        beta: float,
        boundaries: str = ...,
        default_spin: int = ...,
        save_snapshots: bool = ...,
        path_out: Optional[Path] = ...,
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
        """Measure observables and raw acceptance ratios after each thermalization proposal.

        Parameters
        ----------
        N_thermalization : int
            Non-negative number of proposed updates; rejected proposals also count.
        N_resamples : int, optional
            Positive number of bootstrap resamples; use at least two for error estimates.
        observables : list[str]
            Observable names in the order used by all returned arrays.
        seed : int, optional
            Nonzero seeds reproduce a run; zero uses a fresh random seed.
        mu : float
            Hamiltonian parameter (star term).
        h : float
            Hamiltonian parameter (electric field term).
        J : float
            Hamiltonian parameter (plaquette term).
        lmbda : float
            Hamiltonian parameter (gauge field term).
        basis : {'x','z'}, optional
            Spin eigenbasis; selects the C++ backend.
        lattice_type : str
            Geometry: "square", "cubic", "honeycomb", "triangular", or "kagome".
        system_size : int
            Positive linear size; additional restrictions depend on the geometry.
        beta : float
            Positive inverse temperature, equal to the imaginary-time period.
        boundaries : str, optional
            Boundary condition ("periodic", "open").
        default_spin : int, optional
            Default link spin (+1 or -1).
        save_snapshots : bool, optional
            Write GraphML snapshots in addition to returning arrays.
        path_out : pathlib.Path | None, optional
            Existing snapshot directory; None uses the current directory.

        Returns
        -------
        series : numpy.ndarray (complex128), shape (n_obs, N_thermalization)
            Time series per observable during thermalization.
        acc_ratio : numpy.ndarray (float64), shape (N_thermalization,)
            Raw Metropolis ratio (possibly > 1), or zero for an abandoned proposal.
            Acceptance uses min(1, ratio); this array does not record success flags.

        Notes
        -----
        N_resamples is validated but unused by this workflow. Snapshots are recorded
        after proposal 1 and every 10000 proposals thereafter. Real observables have
        zero imaginary part; other entries pack paired real estimators. Returned
        arrays own their data. With no observables, series has shape (0, 0).
        """

    @staticmethod
    def get_sample(
        N_samples: int,
        N_thermalization: int,
        N_between_samples: int,
        *,
        N_resamples: int = ...,
        custom_therm: bool = ...,
        observables: Sequence[str],
        seed: int = ...,
        mu: float,
        h: float,
        h_therm: float = ...,
        J: float,
        lmbda: float,
        lmbda_therm: float = ...,
        basis: str = ...,
        lattice_type: str,
        system_size: int,
        beta: float,
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
        """Thermalize a fresh lattice, then return sampled observables and statistics.

        Parameters
        ----------
        N_samples : int
            Positive number of recorded samples.
        N_thermalization : int
            Non-negative number of initial update proposals.
        N_between_samples : int
            Non-negative number of update proposals before each recorded sample.
        N_resamples : int, optional
            Positive number of bootstrap resamples; use at least two for error estimates.
        custom_therm : bool, optional
            Prepare at (h_therm, lmbda_therm), then ramp h followed by lmbda to the target.
        observables : list[str]
            Observable names in the order used by all returned arrays.
        seed : int, optional
            Nonzero seeds reproduce a run; zero uses a fresh random seed.
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
            Spin eigenbasis; selects the C++ backend.
        lattice_type : str
            Geometry: "square", "cubic", "honeycomb", "triangular", or "kagome".
        system_size : int
            Positive linear size; additional restrictions depend on the geometry.
        beta : float
            Positive inverse temperature, equal to the imaginary-time period.
        boundaries : str, optional
            Boundary condition ("periodic", "open").
        default_spin : int, optional
            Default link spin (+1 or -1).
        save_snapshots : bool, optional
            Write GraphML snapshots in addition to returning arrays.
        path_out : pathlib.Path | None, optional
            Existing snapshot directory; None uses the current directory.

        Returns
        -------
        series : numpy.ndarray (complex128), shape (n_obs, N_samples)
            Time series per observable.
        mean : numpy.ndarray (float64), shape (n_obs,)
            Bias-corrected bootstrap estimate per observable.
        std : numpy.ndarray (float64), shape (n_obs,)
            Bootstrap standard error of the estimate.
        binder : numpy.ndarray (float64), shape (n_obs,)
            Binder ratio; unused estimator categories return zero.
        binder_std : numpy.ndarray (float64), shape (n_obs,)
            Bootstrap error of the Binder ratio.
        tau_int : numpy.ndarray (float64), shape (n_obs,)
            Integrated autocorrelation time in recorded-sample units.

        Notes
        -----
        All series use complex128. Real observables have zero imaginary part; other
        entries pack paired real estimators for Fredenhagen-Marcu or susceptibility.
        The returned arrays own their data. Snapshots include every recorded sample.
        With no requested observables, series has shape (0, 0).
        """

    @staticmethod
    def get_hysteresis(
        N_samples: int,
        N_thermalization: int,
        N_between_samples: int,
        *,
        N_resamples: int = ...,
        observables: Sequence[str],
        seed: int = ...,
        mu: float,
        h_hys: Sequence[float],
        J: float,
        lmbda_hys: Sequence[float],
        basis: str = ...,
        lattice_type: str,
        system_size: int,
        beta: float,
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
        """Traverse paired field schedules once while retaining the lattice state.

        Parameters
        ----------
        N_samples : int
            Positive number of samples per schedule point.
        N_thermalization : int
            Non-negative number of initial update proposals.
        N_between_samples : int
            Non-negative number of update proposals before each recorded sample.
        N_resamples : int, optional
            Positive number of bootstrap resamples; use at least two for error estimates.
        observables : list[str]
            Observable names in the order used by all returned arrays.
        seed : int, optional
            Nonzero seeds reproduce a run; zero uses a fresh random seed.
        mu : float
            Hamiltonian parameter (star term).
        h_hys : list[float]
            Electric field values for each hysteresis step.
        J : float
            Hamiltonian parameter (plaquette term).
        lmbda_hys : list[float]
            Gauge field values for each hysteresis step.
        basis : {'x','z'}, optional
            Spin eigenbasis; selects the C++ backend.
        lattice_type : str
            Geometry: "square", "cubic", "honeycomb", "triangular", or "kagome".
        system_size : int
            Positive linear size; additional restrictions depend on the geometry.
        beta : float
            Positive inverse temperature, equal to the imaginary-time period.
        boundaries : str, optional
            Boundary condition ("periodic", "open").
        default_spin : int, optional
            Default link spin (+1 or -1).
        save_snapshots : bool, optional
            Write GraphML snapshots in addition to returning arrays.
        paths_out : list[pathlib.Path] | None, optional
            Existing directories, one per schedule point when saving snapshots.
            None uses the current directory for every point, overwriting the same file.

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
        The schedules must be nonempty and equally sized. The initial thermalization
        uses h = lmbda = 0, followed by N_thermalization // 4 proposals at each schedule
        point. This call follows the supplied order once; the Python CLI runs both
        forward and reverse branches separately. Series entries pack paired real
        estimators where needed; summary meanings match get_sample(). All arrays own
        their data. With no observables, series3d has shape (n_steps, 0, 0).
        """
