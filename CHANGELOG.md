# ParaToric v1.0.1 Release Notes

Release date: 2026-06-02

v1.0.1 is a minor maintenance release over v1.0 focused on broader
Fredenhagen-Marcu support, faster QMC update paths, leaner result handling, and
refreshed documentation and dependencies.

## Highlights

- Added Wilson and 't Hooft loop construction for honeycomb lattices with
  smooth open boundaries, enabling `fredenhagen_marcu` calculations for that
  geometry.
- Improved Monte Carlo hot paths by reducing unnecessary tuple/star/plaquette
  energy work, avoiding repeated integer-distribution construction, caching edge
  descriptors, and using stack-friendly small buffers for short update paths.
- Reduced avoidable allocation and copying in QMC result construction, Python
  binding conversion, HDF5 output, and hysteresis/sample result assembly.
- Propagated `full_time_series` through sample and hysteresis output handling.
  When full time series output is disabled, sample acceptance-ratio series
  storage is skipped.
- Updated the documentation to the latest arXiv version and clarified "full
  time series" terminology in user-facing text.
- Upgraded the vendored `pybind11` dependency to 3.0.4.

## User-Facing Changes

- The `fredenhagen_marcu` observable now supports the newly implemented open
  honeycomb loop paths.
- The README now uses "Circle-sweep" instead of a LaTeX circle symbol for
  clearer navigation and display.
- Copyright notices were updated to 2026.

## Performance And Internals

- Optimized potential-energy diff calculations by ignoring tuple flips where
  they cannot affect diagonal tuple products.
- Added a direct no-inner-flips edge-energy integral for tuple-flip move
  intervals.
- Added `uniform_index()` for faster unbiased bounded integer sampling.
- Switched several small temporary index, energy, and boolean buffers from
  `std::vector` to `boost::container::small_vector`.
- Moved result buffers instead of copying them in QMC return paths and binding
  conversion paths.
- Added branch-likelihood hints around rare potential-energy resets and
  exceptional update cases.

## Compatibility Notes

- Existing CLI workflows should not require option changes.
- C++ users who explicitly name the return types of low-level lattice/QMC
  energy-difference helpers may need to account for `boost::container::small_vector`
  replacing `std::vector` in those tuple return values. Call sites using `auto`
  should generally be unaffected.
- `--full_time_series 0` continues to write summary statistics while avoiding
  full observable-series output; this release also avoids storing the sample
  acceptance-ratio series in that mode.

## Reference

- Compared with `v1.0`: 13 commits.
- Project-side changes excluding vendored `pybind11`: 22 files changed,
  630 insertions, 272 deletions.
- Source comparison: https://github.com/palmbart/ParaToric/compare/v1.0...v1.0.1

# ParaToric v1.0 Release Notes

Release date: 2026-02-10

v1.0 is the first stable release after v1.0-beta. It expands the observable
set, adds hysteresis workflow support, improves open-boundary percolation
handling, fixes several susceptibility and potential-energy issues, and
strengthens the test suite.

## Highlights

- Added the Python CLI hysteresis sweep workflow and documented
  `etc_hysteresis` usage for forward and backward parameter schedules.
- Added new observables for `largest_plaquette_cluster` and
  `plaquette_percolation_strength`.
- Added support for plaquette percolation probability and strength on open
  boundaries.
- Split susceptibility observables into explicit static and dynamical variants:
  `sigma_x_static_susceptibility`, `sigma_x_dynamical_susceptibility`,
  `sigma_z_static_susceptibility`, and `sigma_z_dynamical_susceptibility`.
- Added integrated autocorrelation-time output and plotting for sweep results.
- Returned/stored acceptance-ratio series for regular sampling runs.
- Added CMake package configuration support via `paratoricConfig.cmake.in`.

## User-Facing Changes

- README examples were refreshed with the v1.0 observable names and larger
  production-oriented sample settings.
- The documentation PDF was updated to a newer arXiv version, and the README now
  includes the arXiv citation link.
- Funding acknowledgements were added to the README.
- Thermalization plots now label the x-axis as "Monte Carlo update" instead of
  "Monte Carlo step".
- Python sweep internals now consistently use `*_lower` and `*_upper` parameter
  names instead of mixed `*_low` and `*_high` names.

## Fixes

- Fixed potential-energy initialization for hysteresis simulations and
  `custom_therm` runs; additional consistency checks now detect cached
  integrated-potential-energy mismatches.
- Fixed high-field/few-flip integrated edge-energy calculations that could
  incorrectly return zero.
- Fixed the dynamical susceptibility factor-of-two issue and adjusted
  susceptibility normalization.
- Fixed combination-update potential-energy cache updates to apply edge
  corrections to the actual affected edges.
- Improved open-style percolation checks so periodic seam links are not treated
  as open-boundary connections.

## Performance And Internals

- Added cached plaquette-to-vertex and star-to-plaquette lookup tables for
  faster star/plaquette combination updates.
- Reworked tuple and edge combination energy-difference paths to reduce
  temporary allocations and avoid unnecessary event-stream work.
- Added `PARATORIC_ENABLE_FAST_MATH`, enabled by default, to compile
  `paratoric_core` with `-ffast-math` on Clang/GNU.
- Improved consistency between star and plaquette update algorithms and cleaned
  up diagnostic output.

## Tests

- Mirrored the `src/` folder structure under `tests/`.
- Added a dedicated `tests/mcmc/test_extended_toric_code_qmc.cpp` suite.
- Added more lattice and QMC test cases across additional lattices and parameter
  regimes.

## Compatibility Notes

- Existing CLI flags remain centered on `*_lower` and `*_upper`, but direct
  Python callers of `JobHandler` sweep methods should use the renamed
  `*_lower`/`*_upper` keyword arguments.
- The old observable names `sigma_x_susceptibility` and
  `sigma_z_susceptibility` were replaced by the explicit static susceptibility
  names. Use the new dynamical names for dynamical susceptibility estimators.
- Builds now use `-ffast-math` by default for Clang/GNU through
  `PARATORIC_ENABLE_FAST_MATH=ON`. Configure with
  `-DPARATORIC_ENABLE_FAST_MATH=OFF` if strict floating-point semantics are
  required.

## Reference

- Compared with `v1.0-beta`: 39 commits.
- Project-side changes excluding vendored dependencies: 15 files changed,
  2150 insertions, 853 deletions.
- Source comparison: https://github.com/palmbart/ParaToric/compare/v1.0-beta...v1.0

# ParaToric v1.0-beta Release Notes

Release date: 2025-10-16

v1.0-beta is the initial public release of ParaToric. It introduced the C++23
continuous-time quantum Monte Carlo implementation for the extended toric code
in parallel fields, together with C/C++/Python interfaces, command line tools,
documentation, and baseline tests.

## Highlights

- Added the core continuous-time QMC implementation for the extended/perturbed
  toric code at finite temperature in the `x` and `z` bases.
- Added lattice construction and update support for square, cubic, triangular,
  honeycomb, and kagome-style geometries, with periodic/open boundary handling
  where implemented.
- Added public C++ APIs via `paratoric::ExtendedToricCode` for thermalization,
  sampling, and hysteresis simulations.
- Added a C API wrapper with explicit status codes, result structs, and result
  cleanup support.
- Added optional Python bindings through vendored `pybind11`, with typed package
  metadata and NumPy-based return values.
- Added C++ and Python command line workflows for sample, thermalization,
  temperature sweeps, field sweeps, circle sweeps, and hysteresis runs.
- Added HDF5 output for observables, summary statistics, Binder ratios,
  autocorrelation estimates, acceptance ratios, and optional full time series.

## Observable Set

- Included energy observables: `energy`, `energy_h`, `energy_mu`, `energy_J`,
  and `energy_lmbda`.
- Included local and topological observables such as `sigma_x`, `sigma_z`,
  `star_x`, `plaquette_z`, `delta`, `string_number`, `anyon_count`,
  `anyon_density`, `staggered_imaginary_times`, and `fredenhagen_marcu`.
- Included percolation observables for string/plaquette/cube diagnostics,
  including `percolation_probability`, `percolation_strength`,
  `plaquette_percolation_probability`, and `cube_percolation_probability`.
- Included the beta-era susceptibility names `sigma_x_susceptibility` and
  `sigma_z_susceptibility`.

## Tooling And Packaging

- Added a CMake build for `paratoric_core`, the optional C++ CLI, tests, and
  optional Python bindings.
- Added support for native CPU optimization, AVX2 flags, and optional MPI
  linking.
- Added Python packaging metadata for the `paratoric` package.
- Added a SLURM job-script creation notebook under `scripts/`.
- Added project documentation PDF and logo assets.
- Added Boost unit tests for lattice construction and extended toric code QMC
  behavior.

## Known Beta Limitations

- The `v1.0-beta` tag was the first public snapshot, so there is no prior
  release comparison.
- The CMake package config template referenced by the build was missing in this
  beta tag; it was added in `v1.0`.
- Hysteresis support was present in the APIs and Python CLI, but README
  documentation for the workflow was added later in `v1.0`.
- Dynamical susceptibility observables were partially implemented internally but
  not exposed in the beta observable map.
- Open-boundary honeycomb Wilson loops and honeycomb/kagome 't Hooft loop
  support were not yet available for all Fredenhagen-Marcu use cases.
- Some percolation and susceptibility estimators were refined in the subsequent
  `v1.0` release.

## Reference

- Initial public release commit: `5fb6b7d`.
- Initial snapshot size: 392 files changed, 99195 insertions.
- Source tag: https://github.com/palmbart/ParaToric/releases/tag/v1.0-beta
