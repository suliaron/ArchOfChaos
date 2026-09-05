# Changelog

All notable changes to **Arch of Chaos** are documented in this file.

The project follows semantic-style version numbering where practical. Git tags and GitHub Releases should be used to preserve the exact source state of each published version.

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
