# Changelog

All notable changes to **Arch of Chaos** are documented in this file.

The project follows semantic-style version numbering where practical. Git tags and GitHub Releases should be used to preserve the exact source state of each published version.

---

## [1.7.0]

### Added

- Comprehensive reproducibility header for all numerical output files.
- Program and run metadata written automatically to every output file, including:
  - program name and version;
  - author and affiliation;
  - run timestamp;
  - build timestamp;
  - build configuration;
  - compiler;
  - C++ standard;
  - operating-system/platform information;
  - host/computer name;
  - numerical integration method;
  - run mode;
  - chaos indicator;
  - mathematical formalism.
- Exact commented copy of the complete input file embedded in every output file.
- Machine-readable description of the numerical output structure.
- Common `DATA` section separating metadata from numerical results.
- Centralized numerical output formatting and table-header generation in the I/O module.
- New logarithmic output scheduler for `INDICATOR` mode.
- New `LogOutputSchedule` class in `time_utils.h` / `time_utils.cpp`.
- New INDICATOR input parameter:

  ```text
  output_first
  ```

  defining the first elapsed physical output time.

- New INDICATOR input parameter:

  ```text
  output_points_per_decade
  ```

  defining the number of output points in each time decade.

- Default logarithmic INDICATOR output configuration:

  ```text
  output_first = 1.0
  output_points_per_decade = 9
  ```

- Output-schedule metadata in the reproducibility header:
  - linear output interval for `ORBIT`;
  - logarithmic first-output time and points per decade for `INDICATOR`;
  - final-value-only output for `GRID`.

### Changed

- `INDICATOR` mode now uses logarithmically distributed output times instead of a fixed linear `output_dt`.
- `output_dt` is now used exclusively by `ORBIT` mode.
- INDICATOR output times are distributed linearly inside each decade and the decade scale increases by a factor of ten.

  For example,

  ```text
  output_first = 1.0e-4
  output_points_per_decade = 9
  ```

  produces

  ```text
  1e-4, 2e-4, ..., 9e-4,
  1e-3, 2e-3, ..., 9e-3,
  1e-2, ...
  ```

- Arbitrary positive first-output times are supported, including sub-day values.
- The final integration point is always written even if it does not coincide with a scheduled output time.
- A final integration point that already coincides with a scheduled output time is written only once.
- The initial FLI value at `t = 0` is no longer written in INDICATOR mode, allowing the result to be plotted directly on logarithmic time axes.
- LCI evaluation in INDICATOR mode now uses physical elapsed time, preserving the unit `1/day`.
- Automatic INDICATOR output-file names now describe the logarithmic schedule.

  Example:

  ```text
  FLI_NEWTONIAN_T-0.01_log-0.0001-N9_...
  ```

  instead of the obsolete fixed-interval form:

  ```text
  FLI_NEWTONIAN_T-0.01_dt-0_...
  ```

- Automatic ORBIT filenames retain the existing `_dt-<value>` convention.
- Output-table formatting has been removed from individual run functions and centralized in the I/O module.
- Output files now use common scientific formatting with explicit signs and consistent field widths.
- Planar inclination validation now explicitly requires

  ```text
  0 <= i < PLANAR_EPS
  ```

  instead of using an absolute-value test.
- The same inclination restriction is applied consistently to:
  - fixed ORBIT initial conditions;
  - fixed INDICATOR initial conditions;
  - fixed GRID inclinations;
  - grid-controlled inclination ranges.

### Validation and regression tests

- Verified `LogOutputSchedule` with:

  ```text
  output_first = 1.0e-4
  output_points_per_decade = 9
  ```

- Verified a non-default schedule with:

  ```text
  output_first = 1.0e-4
  output_points_per_decade = 5
  ```

- Verified logarithmic INDICATOR output across multiple time decades.
- Verified that a non-scheduled final integration time is always written.
- Verified that a scheduled final integration time is not duplicated.
- Verified automatic INDICATOR filename generation with logarithmic scheduling parameters.
- Verified reproducibility-header output including the complete input-file copy and numerical output description.
- Successfully completed a 100-year FLI `(a,e)` GRID calculation with 201201 initial conditions.
- Verified the version update and command-line version output for Arch of Chaos 1.7.0.

---

## [1.6.0]

### Added

- Generalized GRID mode from the fixed `(a,e)` plane to arbitrary multidimensional orbital-element grids.
- New generic grid input syntax:

  ```text
  grid = <orbital element> <minimum> <maximum> <number of intervals>
  ```

- Support for the grid elements:

  ```text
  a, e, i, omega, Omega, tau, M
  ```

- `OrbitalElement` and `GridAxis` abstractions.
- Generic multidimensional `GridIterator`.
- Element-based grid lookup and `hasAxis()` support.
- `GridIterator::apply()` for applying the current grid point to orbital elements.
- Periodic endpoint handling for full angular ranges such as `0 ... 360 deg`.
- Generalized GRID output headers generated from the input grid axes.
- LCI output unit shown explicitly as `LCI [1/day]`.
- Generalized GRID integration using the new grid iterator.
- Grid-controlled mean anomaly handling with conversion from `M` to `tau` using the current semimajor axis.
- Per-grid-point recalculation of physical integration duration when `nPeriods` is used.
- Automatic output-file naming for `GRID`, `ORBIT`, and `INDICATOR` modes.
- Common time-series filename generation for `ORBIT` and `INDICATOR`.
- Automatic GRID filenames containing grid-axis limits and interval counts.
- Verbose reporting of the actual output path.
- Additional read-only `InitData` accessors required by generalized execution and automatic filename generation.

### Changed

- GRID mode is no longer restricted to the fixed `(a,e)` parameter plane.
- Grid axes now preserve their input-file order; the first axis varies fastest.
- Fixed and grid-controlled orbital elements are validated as mutually exclusive sources.
- `tau` and `M` are treated consistently as alternative representations of orbital phase.
- Full-period angular grids omit the duplicated upper endpoint.
- `T` is treated as a physical integration duration in days.
- LCI in GRID mode is normalized using physical elapsed time, preserving the unit `1/day`.
- GRID verbose output now lists each grid axis together with its limits, interval count, and unit.
- Output files can now be generated without explicitly supplying `-o`.
- Explicit `-o` output still has priority over automatic naming.

### Validation and regression tests

- Verified multidimensional grid traversal.
- Verified reversed grid-axis order.
- Verified three-dimensional `(a,e,M)` traversal.
- Verified periodic `M = 0 ... 360 deg` endpoint handling.
- Verified partial `M = 0 ... 180 deg` endpoint handling.
- Verified total grid-point counting and multidimensional index carry.
- Verified `GridIterator::apply()` and degree-to-radian conversion.
- Verified preservation of non-grid-controlled orbital elements.
- Verified small `(a,e)` FLI grids with the generalized GRID implementation.
- Verified small `(a,e)` LCI grids with the generalized GRID implementation.
- Verified automatic output naming for GRID mode.
- Verified automatic output naming for ORBIT mode.
- Verified automatic output naming for INDICATOR mode.
- Verified explicit `-o` override behavior.

---

## [1.5.0]

Previous stable development milestone before the generalized orbital-element grid and automatic output-file naming work introduced in 1.6.0.

For the exact historical source state, use the corresponding Git tag/release once published.
