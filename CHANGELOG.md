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
