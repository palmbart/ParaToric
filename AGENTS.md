# ParaToric: development and integration guide

ParaToric implements continuous-time quantum Monte Carlo for the extended toric
code in parallel fields, in the spin `x` or `z` basis. The C++23 core supplies C++,
C, and Python APIs; a native CLI writes HDF5, and a separate Python CLI runs
parameter sweeps and plots results. License: MIT. Publication citations and
longer usage examples are in `README.md`; the full manual is
`doc/Documentation.pdf`.

Use this file as the starting map. Read the relevant source/header rather than
loading the entire repository or PDF. CMake files and public headers are the
authority for current build and API behavior; update this guide when those change.

## Source map

| Location | Responsibility |
| --- | --- |
| `include/paratoric/mcmc/extended_toric_code.hpp` | Public C++ facade: `paratoric::ExtendedToricCode`. |
| `include/paratoric/types/types.hpp` | Public `Config`, its four specification structs, and `Result`. |
| `include/paratoric/mcmc/extended_toric_code_c.h` | Public C ABI, ownership, and status codes. |
| `src/mcmc/extended_toric_code.cpp`, `extended_toric_code_c.cpp` | Runtime basis dispatch and C adapter. |
| `src/mcmc/extended_toric_code_qmc.hpp` | Template backend `ExtendedToricCodeQMC<'x'/'z'>`, updates, estimators, `obs_vec` registry. |
| `src/mcmc/input_validation.hpp` | Shared configuration validation; geometry/coupling checks also live in lattice/backend code. |
| `src/lattice/lattice.{hpp,cpp}`, `time_search.hpp` | Geometry, event histories, energy caches, percolation/loops, GraphML snapshots, sorted time searches. |
| `src/rng/`, `src/statistics/` | RNG, bootstrap estimators, autocorrelation. |
| `src/cli/paratoric.cpp`, `src/io/` | Boost.Program_options native CLI and HDF5 serialization. |
| `python/bindings/`, `python/paratoric/` | pybind11 extension, package imports, `_paratoric.pyi` signatures and array contracts. |
| `python/cli/paratoric.py`, `job_handler.py` | Sweep arguments, multiprocessing, native subprocesses, HDF5 aggregation, plots. |
| `tests/` | Boost.Test executables registered with CTest. |
| `CMakeLists.txt`, `cmake/`, `python/pyproject.toml` | Targets/options, exported CMake package, Python metadata. |
| `external/pybind11/`, `scripts/job_script_creator.ipynb` | Optional binding submodule; SLURM job generation notebook. |

Only the three headers under `include/paratoric/` are installed public APIs.
Headers under `src/` are private implementation details.

## Build and test

Run commands from the checkout root. Requirements: CMake >= 3.23 and a C++23
compiler **and standard library** (including `std::format`). README documents
GCC 15 / Clang 20, Boost >= 1.87, and HDF5 >= 1.14.3; CMake does not enforce the
Boost/HDF5 version minima.

- Core: Boost `log_setup`, `log` (CMake config packages).
- Native CLI: additionally Boost `program_options` and HDF5 C and C++ libraries.
- Tests: additionally Boost `unit_test_framework`.
- Python: >= 3.10, NumPy >= 1.24, matplotlib >= 3.6, h5py >= 3.7. Bindings also
  need Python development headers and the `external/pybind11` submodule.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/install"
cmake --build build -j4
ctest --test-dir build -j4 --output-on-failure
cmake --install build
```

Set an explicit writable install prefix: current CMake uses its normal default
(typically `/usr/local`), despite the README's older checkout-prefix claim.
Installation exports headers, `libparatoric_core.a`, and the `paratoric` CMake
package under the prefix. A top-level CLI install **also copies the executable
to the checkout's `bin/paratoric`**, which the Python CLI uses.

| CMake option | Default / effect |
| --- | --- |
| `PARATORIC_BUILD_CLI`, `PARATORIC_BUILD_TESTS` | ON for top-level builds, OFF when embedded with `add_subdirectory`. |
| `PARATORIC_BUILD_PYBIND` | OFF; builds `_paratoric` when ON. |
| `PARATORIC_ENABLE_NATIVE_OPT`, `PARATORIC_ENABLE_AVX2` | OFF; opt into CPU-specific flags for suitable machines. |
| `PARATORIC_ENABLE_FAST_MATH` | ON; set OFF when checking strict floating-point behavior. |
| `PARATORIC_LINK_MPI` | OFF; attempts optional MPI linkage, continues without MPI if unavailable. Does not parallelize the sampler. |
| `PARATORIC_EXPORT_COMPILE_COMMANDS` | ON; generates compile commands with supported generators. |

For core-only builds, add `-DPARATORIC_BUILD_CLI=OFF -DPARATORIC_BUILD_TESTS=OFF`;
HDF5, Python, and pybind11 are then unnecessary. The core is explicitly STATIC;
`BUILD_SHARED_LIBS` does not change it. For dependency discovery use
`-DCMAKE_PREFIX_PATH=/dependency/prefix` (or `Boost_DIR` / `HDF5_DIR`). Select a
compiler with `CC`/`CXX` or `CMAKE_CXX_COMPILER` in a fresh build directory.

For development use a separate Debug build, e.g. `-B build/debug
-DCMAKE_BUILD_TYPE=Debug`. CTest names are `test_lattice`,
`test_extended_toric_code_qmc`, and `test_input_validation`; select a relevant
test with `ctest --test-dir build/debug -R test_input_validation --output-on-failure`.

## C++ and C integration

Installed package consumer (`main.cpp` below):

```cmake
cmake_minimum_required(VERSION 3.23)
project(example LANGUAGES C CXX)
find_package(Boost REQUIRED COMPONENTS log_setup log CONFIG)
find_package(paratoric CONFIG REQUIRED)
add_executable(example main.cpp)
target_link_libraries(example PRIVATE paratoric::core)
```

Configure the consumer with `-DCMAKE_PREFIX_PATH=/path/to/install/prefix`, or use
`-Dparatoric_DIR=/path/to/paratoric/build/cmake` to consume an already-built tree.
The package config currently **does not call `find_dependency`**: consumers must
find Boost themselves, as above. If the core was built with MPI linkage, also
call `find_package(MPI REQUIRED COMPONENTS C CXX)`.

For vendoring, replace both `find_package` lines with
`add_subdirectory(path/to/paratoric paratoric-build)`, then link the same target.
The target propagates public includes, Boost linkage, and C++23 requirements.
Set `CMAKE_POSITION_INDEPENDENT_CODE=ON` when embedding the static core in a
shared library (the Python binding build enables PIC automatically).

Minimal C++ smoke example (small counts exercise the API, not convergence):

```cpp
#include <paratoric/mcmc/extended_toric_code.hpp>

int main() {
    paratoric::Config cfg{};
    cfg.lat_spec.system_size = 4;
    cfg.lat_spec.beta = 1.0;
    cfg.sim_spec.N_thermalization = 100;
    cfg.sim_spec.N_samples = 20;
    cfg.sim_spec.N_between_samples = 10;
    cfg.sim_spec.N_resamples = 20;
    cfg.sim_spec.seed = 123;
    cfg.sim_spec.observables = {"energy"};
    auto result = paratoric::ExtendedToricCode::get_sample(cfg);
    return result.mean.size() == 1 ? 0 : 1;
}
```

All three static methods (`get_sample`, `get_thermalization`, `get_hysteresis`)
take `const Config&` and return `Result`. Each call creates a fresh backend;
state persists across parameter points only within one hysteresis call.
Invalid inputs throw `std::invalid_argument`; runtime failures throw exceptions.

For C, include `<paratoric/mcmc/extended_toric_code_c.h>` and link the same core
with a C++ linker (CMake: enable CXX and set the C executable's `LINKER_LANGUAGE`
to `CXX`). Fill `ptc_config_t` explicitly; zero initialization alone is not a
valid simulation. Use `ptc_create`, `ptc_get_sample` / `ptc_get_thermalization` /
`ptc_get_hysteresis`, and check `ptc_status_t` plus `ptc_last_error()`. Inputs are
borrowed during the synchronous call. Zero-initialize `ptc_result_t`, release
its allocations with `ptc_result_destroy` before reuse, and release the handle
with `ptc_destroy`.

## Python API and CLIs

In an activated virtual environment, build the Python API separately:

```sh
git submodule update --init --recursive external/pybind11
python -m pip install -e ./python
cmake -S . -B build/python -DPARATORIC_BUILD_PYBIND=ON \
  -DPARATORIC_BUILD_CLI=OFF -DPARATORIC_BUILD_TESTS=OFF \
  -DPython_EXECUTABLE="$VIRTUAL_ENV/bin/python" \
  -DPython3_EXECUTABLE="$VIRTUAL_ENV/bin/python"
cmake --build build/python -j4
python -c "from paratoric import extended_toric_code; print(extended_toric_code.get_sample.__doc__)"
```

`pip install -e ./python` installs dependencies/metadata; it **does not compile
the extension**. CMake writes `_paratoric*.so` / `.pyd` into
`python/paratoric/`. Import `from paratoric import extended_toric_code` and call
the three `get_*` functions; see `_paratoric.pyi` for required keyword arguments.
`get_sample` returns `(series, mean, std, binder, binder_std, tau_int)` as NumPy
arrays. `series` has shape `(n_obs, N_samples)` and dtype `complex128`; summary
arrays are `float64`. Hysteresis adds a leading schedule-point dimension.
Thermalization returns `(series, acc_ratio)`. Invalid inputs raise `ValueError`.

The two CLIs have different workflows:

- `build/paratoric --help`: native `etc_sample`, `etc_thermalization`,
  `etc_hysteresis`; takes `--beta`. Writes per-run `obs.h5` and optional GraphML
  snapshots. Hysteresis requires one `--folder_names` entry per schedule point.
- `python python/cli/paratoric.py --help`: `etc_T_sweep`, `etc_h_sweep`,
  `etc_lmbda_sweep`, `etc_circle_sweep`, `etc_hysteresis`, `etc_thermalization`;
  uses temperature `T=1/beta`, multiprocessing and the installed
  **checkout `bin/paratoric`**, not the extension or a PATH lookup. Needs NumPy,
  matplotlib, h5py but no binding build. Produces `parameters.txt`, plots,
  aggregate `simulation_data.h5`, and raw per-job `obs.h5` files.

`python -m paratoric` launches the Python CLI in an editable/source checkout;
it depends on sibling `python/cli/` files, which normal package discovery does
not include. See README for full sweep commands.

## Simulation contracts and development checks

- Couplings: `mu` = star, `J` = plaquette, `h` = x field, `lmbda` = z field.
  Off-diagonal couplings must be nonnegative: `J,lmbda` in x basis; `mu,h` in z.
- Geometries: `square`, `cubic`, `triangular`, `honeycomb`, `kagome`; boundaries
  `periodic`/`open`. Periodic triangular/honeycomb/kagome require even size.
  Observable support depends on geometry; consult lattice implementations.
- `beta > 0`, size > 0, spin +/-1; sample/resample counts > 0, thermalization and
  spacing counts >= 0. Counts refer to update proposals, including rejections.
  Use at least two bootstrap resamples for error estimates. Nonzero seeds
  reproduce a run; public calls with seed 0 initialize a fresh random RNG.
- Observable order is preserved. C++ series entries are real/complex variants;
  complex values can pack paired real estimators for Fredenhagen-Marcu and
  susceptibilities. Preserve their statistics categories. `acc_ratio` stores
  raw Metropolis ratios (possibly > 1), not accepted/rejected flags.
- Hysteresis schedules must be nonempty and equally sized. The API traverses
  them once; the Python CLI runs forward and reverse branches separately.
- Library APIs return data in memory; HDF5 belongs to the CLI layer.
  `full_time_series` controls CLI serialization, not API return values.
  Create snapshot directories before API calls that enable snapshot output.
- For backend changes preserve x/z duality, sorted imaginary-time histories,
  and consistency of accepted moves with energy caches. Rebuild bare caches
  when changing couplings; zero-coupling caches may be stale.
- Add observables in backend `obs_vec` and Python `JobHandler.obs_dict`; update
  statistics/I/O handling if introducing a category. For geometry changes,
  update validation and check loops/percolation as well as lattice construction.
- Keep C++, C, pybind11, and Python stubs consistent when changing public
  contracts. Build the affected interfaces and run relevant CTest cases;
  exercise CLI/Python paths separately when changing them. QMC smoke tests do
  not establish physical convergence. For documentation-only changes, check
  paths and examples without requiring the entire simulation suite.
- Keep builds, binaries, simulation output, and virtual environments out of
  source commits. Preserve unrelated working-tree edits and avoid modifying
  the pybind11 submodule unless the task specifically requires it.
