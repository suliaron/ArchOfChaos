# Arch of Chaos

**Arch of Chaos** is a C++ numerical tool for studying orbital dynamics and chaos in the **planar circular restricted three-body problem (CRTBP)**.

The program integrates trajectories of a massless third body in the gravitational field of two primaries moving on circular orbits. It supports both **Newtonian** and **Hamiltonian** formulations, computes chaos indicators from the variational equations, and can evaluate those indicators over arbitrary multidimensional grids of orbital elements.


> **Project status:** active research/development software.  
> The current implementation has been tested against selected regression cases, but it should still be treated as research code rather than a general-purpose production library.

---

## Contents

- [Overview](#overview)
- [Main features](#main-features)
- [Scientific model](#scientific-model)
- [Units and conventions](#units-and-conventions)
- [Numerical integration](#numerical-integration)
- [Chaos indicators](#chaos-indicators)
- [Run modes](#run-modes)
- [Building](#building)
- [Command-line interface](#command-line-interface)
- [Input-file format](#input-file-format)
- [Integration duration](#integration-duration)
- [Orbital phase](#orbital-phase)
- [GRID mode](#grid-mode)
- [Output files](#output-files)
- [Automatic output-file naming](#automatic-output-file-naming)
- [Examples](#examples)
- [Plotting GRID results](#plotting-grid-results)
- [Validation and regression tests](#validation-and-regression-tests)
- [Project structure](#project-structure)
- [Known limitations](#known-limitations)
- [Suggested repository additions](#suggested-repository-additions)

---

# Overview

Arch of Chaos is designed for numerical experiments in celestial mechanics where one needs to:

- integrate a single CRTBP trajectory;
- compare Newtonian and Hamiltonian formulations;
- integrate the variational equations together with the orbit;
- compute finite-time chaos indicators;
- generate high-resolution dynamical maps;
- scan arbitrary combinations of osculating orbital elements;
- reproduce numerical experiments with explicitly stored physical and numerical parameters.

The program uses physical orbital elements as input, constructs the corresponding heliocentric inertial Cartesian state, transforms it to the dimensionless rotating CRTBP frame, and then performs the numerical integration.

A typical workflow is

```text
input file
    |
    v
InitData parser + validation
    |
    v
orbital elements
    |
    v
heliocentric inertial Cartesian state
    |
    v
dimensionless rotating CRTBP state
    |
    v
Newtonian or Hamiltonian equations
    |
    v
RKF54 integration
    |
    +--> ORBIT time series
    |
    +--> FLI / LCI time series
    |
    +--> multidimensional GRID map
```

---

# Main features

## Dynamical model

- Planar circular restricted three-body problem.
- Two massive primaries, `P1` and `P2`.
- Massless third body `P3`.
- Physical input in astronomical units, solar masses, and days.
- Internal transformation to a normalized rotating CRTBP system.

## Mathematical formulations

- `NEWTONIAN`
- `HAMILTONIAN`

The two formulations have been explicitly compared in development tests.

## Chaos indicators

Currently implemented:

- **FLI** — Fast Lyapunov Indicator
- **LCI** — Lyapunov Characteristic Indicator

Recognized but not yet implemented:

- **RLI** — Relative Lyapunov Indicator

## Computation modes

- `ORBIT`
- `INDICATOR`
- `GRID`

## Generalized orbital-element grids

GRID mode supports arbitrary combinations of

```text
a
e
i
omega
Omega
tau
M
```

using the syntax

```text
grid = <orbital element> <minimum> <maximum> <number of intervals>
```

The order in which grid axes are written in the input file determines the iteration order. The **first grid axis varies fastest**.

## Integration-duration input

The integration duration can be specified either as

```text
T = <physical duration in days>
```

or

```text
nPeriods = <number of orbital periods of P3>
```

These are alternative representations of the integration duration.

## Orbital-phase input

The initial orbital phase can be specified either as

```text
tau = <time of pericenter passage in days>
```

or

```text
M = <mean anomaly in degrees>
```

These are alternative representations of the same orbital-phase information.

## Automatic output naming

If no `-o` command-line option is supplied, the program automatically constructs a descriptive output file name from the model parameters.

Explicit `-o` output always has priority over the generated name.

---

# Scientific model

Let the two primary masses be

\[
m_1,\qquad m_2,
\]

and let their constant separation be

\[
a_2.
\]

The CRTBP mass parameter is

\[
\mu = \frac{m_2}{m_1+m_2}.
\]

The physical gravitational parameter of the relative `P1-P2` orbit is

\[
\mu_{12}=k^2(m_1+m_2),
\]

where `k` is the Gaussian gravitational constant in the unit system used by the program.

The mean motion of the primary system is

\[
n = \sqrt{\frac{\mu_{12}}{a_2^3}}.
\]

For the heliocentric `P1-P3` Keplerian orbit,

\[
\mu_{13}=k^2 m_1.
\]

The input osculating orbital elements of `P3` are

\[
(a,e,i,\omega,\Omega,\tau)
\]

or equivalently

\[
(a,e,i,\omega,\Omega,M).
\]

The orbital elements are first converted to a heliocentric inertial Cartesian state and then transformed to the rotating normalized CRTBP coordinates.

---

# Units and conventions

The input file uses the following physical units.

| Quantity | Input unit |
|---|---|
| `m1`, `m2` | solar mass |
| `a2` | AU |
| semimajor axis `a` | AU |
| eccentricity `e` | dimensionless |
| `i` | degree |
| `omega` | degree |
| `Omega` | degree |
| `M` | degree |
| `tau` | day |
| `t0` | day |
| `T` | day |
| `output_dt` | day |
| `relTol`, `absTol` | dimensionless numerical tolerances |

Angular input values are converted internally to radians.

The numerical CRTBP integration is performed in **dimensionless CRTBP time**. Physical times are converted to and from the normalized system as needed.

For LCI output, the final normalization uses **physical elapsed time in days**, therefore the LCI unit is

```text
1/day
```

---

# Numerical integration

Arch of Chaos currently uses an adaptive **Runge-Kutta-Fehlberg 5(4)** integrator (`RKF54`).

The integration machinery includes:

- adaptive step-size control;
- relative tolerance `relTol`;
- absolute tolerance `absTol`;
- exact limiting of the final integration step;
- periodic finite-state checks;
- integration of the variational equations when a chaos indicator is requested.

A typical tolerance choice is

```text
relTol = 1.0e-6
absTol = 1.0e-10
```

These values are examples only; appropriate tolerances depend on the intended experiment.

---

# Chaos indicators

## FLI

The Fast Lyapunov Indicator is evaluated from the deviation-vector norm.

The implemented form is a running maximum,

\[
\mathrm{FLI}(t)
=
\max_{\tau\le t}
\log \|\delta(\tau)\|.
\]

The deviation vector is integrated together with the orbit.

In GRID mode, one final FLI value is written for every grid point.

---

## LCI

The Lyapunov Characteristic Indicator is evaluated as

\[
\mathrm{LCI}(t)
=
\frac{1}{t-t_0}
\ln\left(
\frac{\|\delta(t)\|}
     {\|\delta(t_0)\|}
\right).
\]

For the current physical output convention,

\[
t-t_0
\]

is expressed in **days**, and therefore

\[
[\mathrm{LCI}] = \mathrm{day}^{-1}.
\]

In GRID mode, LCI is evaluated at the final integration time for every grid point.

---

## RLI

`RLI` is recognized by the input/parser infrastructure but is **not yet implemented**.

Attempting to run an unsupported RLI calculation results in an error.

---

# Run modes

## ORBIT

`ORBIT` integrates one initial condition and writes the time evolution of the orbit.

Use

```text
mode = ORBIT
indicator = NONE
```

The initial osculating orbital elements must be supplied as fixed values.

The output sampling interval is controlled by

```text
output_dt = ...
```

---

## INDICATOR

`INDICATOR` integrates one initial condition together with the variational equations and writes the selected chaos indicator as a function of time.

Example:

```text
mode = INDICATOR
indicator = FLI
formalism = NEWTONIAN
```

or

```text
mode = INDICATOR
indicator = LCI
formalism = HAMILTONIAN
```

`output_dt` determines the output cadence.

---

## GRID

`GRID` evaluates one final chaos-indicator value at every point of a multidimensional orbital-element grid.

Example:

```text
mode = GRID
indicator = LCI
formalism = NEWTONIAN
```

The output columns are generated automatically from the grid axes.

---

# Building

The project is currently developed and tested on **Windows with Microsoft Visual Studio**.

A typical tested configuration is

```text
Configuration: Release
Platform:      x64
```

The executable is typically produced as

```text
archofchaos.exe
```

For example, a development build may appear under a path such as

```text
x64\Release\archofchaos.exe
```

The exact C++ language standard and compiler options are controlled by the Visual Studio project settings.

At present this README does not define a separate CMake/Linux build workflow.

---

# Command-line interface

Basic usage:

```text
archofchaos -i <input file>
```

Explicit output file:

```text
archofchaos -i <input file> -o <output file>
```

Verbose mode:

```text
archofchaos -i <input file> --verbose
```

Common options:

```text
-i <file>       input initialization file
-o <file>       explicit output file
-h, --help      display help
-v, --version   display program version
--verbose       print initialization and run information
```

If `-o` is omitted, Arch of Chaos generates an output file name automatically from the input parameters.

With `--verbose`, the program also reports the actual output path.

Example in PowerShell:

```powershell
& "D:\path\to\archofchaos.exe" `
    -i .\input.txt `
    --verbose
```

---

# Input-file format

The input is a plain-text file containing

```text
key = value
```

assignments.

Comments start with `#`.

Example:

```text
# Example Arch of Chaos input

mode = ORBIT
indicator = NONE
formalism = NEWTONIAN

m1 = 1.0
m2 = 0.001001001001001
a2 = 5.2026

t0 = 0.0
nPeriods = 10

relTol = 1.0e-6
absTol = 1.0e-10

output_dt = 10.0

a = 5.2
e = 0.1
i = 0.0
omega = 60.0
Omega = 0.0
M = 30.0
```

---

# Integration duration

Exactly one integration-duration representation should be supplied.

## Physical duration

```text
T = 1000.0
```

means

```text
integrate for 1000 days
```

`T` is a **duration**, not an absolute final epoch.

If

```text
t0 = 500.0
T  = 1000.0
```

then the corresponding physical final epoch is

\[
t_\mathrm{final}=1500\ \mathrm{day}.
\]

---

## Number of P3 orbital periods

```text
nPeriods = 100
```

means that the orbit is integrated for 100 Keplerian orbital periods of `P3`.

For the heliocentric `P1-P3` orbit,

\[
n_3 = \sqrt{\frac{\mu_{13}}{a^3}},
\]

and

\[
P_3 = \frac{2\pi}{n_3}.
\]

Therefore

\[
T = n_\mathrm{Periods} P_3.
\]

### Important GRID behavior

If `a` is a grid axis and `nPeriods` is used, the physical integration duration is recalculated **independently for every grid point** using the current semimajor axis.

Thus all grid points are integrated for the same number of their own orbital periods, not necessarily for the same number of physical days.

---

# Orbital phase

The orbital phase can be supplied by either `tau` or `M`.

## Time of pericenter passage

```text
tau = 0.0
```

Unit:

```text
day
```

## Mean anomaly

```text
M = 60.0
```

Unit:

```text
degree
```

The relation used is

\[
M(t_0)=n_3(t_0-\tau),
\]

so that

\[
\tau=t_0-\frac{M(t_0)}{n_3}.
\]

In GRID mode, if `M` is grid-controlled, `tau` is recalculated using the **current grid value of `a`**.

`tau` and `M` are mutually exclusive phase representations.

---

# GRID mode

## Generic grid syntax

A grid axis is defined as

```text
grid = <orbital element> <minimum> <maximum> <number of intervals>
```

Example:

```text
grid = a 4.8 5.8 250
grid = M 30.0 330.0 30
```

Supported names:

```text
a
e
i
omega
Omega
tau
M
```

The number given at the end is the **number of intervals**, not the number of points.

For an ordinary non-periodic grid,

```text
grid = e 0.0 0.2 2
```

contains

```text
0.0
0.1
0.2
```

and therefore has 3 grid points.

---

## Grid iteration order

The first axis listed in the input varies fastest.

Example:

```text
grid = a 3.0 4.0 1
grid = e 0.0 0.2 2
```

is traversed as

```text
a=3  e=0.0
a=4  e=0.0

a=3  e=0.1
a=4  e=0.1

a=3  e=0.2
a=4  e=0.2
```

If the order is reversed,

```text
grid = e 0.0 0.2 2
grid = a 3.0 4.0 1
```

then `e` varies fastest.

---

## Periodic angular endpoints

For a complete angular period, the duplicated upper endpoint is omitted.

Example:

```text
grid = M 0.0 360.0 4
```

produces

```text
0
90
180
270
```

and not

```text
0
90
180
270
360
```

because `0 deg` and `360 deg` represent the same direction.

For a partial interval,

```text
grid = M 0.0 180.0 18
```

both endpoints are retained:

```text
0, 10, 20, ..., 170, 180
```

---

## Fixed versus grid-controlled orbital elements

In GRID mode each orbital quantity must have a unique source.

An element may be

```text
fixed
```

or

```text
grid-controlled
```

but not both.

Invalid:

```text
a = 5.2
grid = a 4.8 5.8 100
```

The program rejects this because `a` is specified twice.

Likewise, an orbital element required by the initial condition cannot simply be omitted unless it is supplied as a grid axis.

---

## Phase validation in GRID mode

The following rules apply:

- `tau` and `M` cannot both be grid axes.
- A fixed phase and a grid-controlled phase cannot be specified simultaneously.
- If neither `tau` nor `M` is a grid axis, one of them must be supplied as a fixed input value.
- `tau` and `M` remain alternative representations of orbital phase.

---

## Planar-model restriction

Although the generic grid infrastructure recognizes `i`, the current dynamical model is planar.

Therefore the planar CRTBP validation requires

```text
i = 0
```

within numerical tolerance.

A genuinely three-dimensional spatial CRTBP is not currently implemented.

---

# Output files

Every output file starts with a program/version comment, for example

```text
# Arch of Chaos 1.6.0
```

---

## ORBIT output

ORBIT mode writes a time series of the dynamical state.

The exact columns depend on the current output implementation and mathematical formalism, but the physical output times are expressed in days.

---

## INDICATOR output

FLI example:

```text
# Arch of Chaos 1.6.0
t [day]           FLI
...
```

LCI example:

```text
# Arch of Chaos 1.6.0
t [day]           LCI [1/day]
...
```

---

## GRID output

GRID columns follow the input grid-axis order.

Example for an `(a,M)` map:

```text
# Arch of Chaos 1.6.0
a [AU]            M [deg]           LCI [1/day]
4.8000000000e+00  3.0000000000e+01  2.0686669529e-03
4.8040000000e+00  3.0000000000e+01  2.0669752898e-03
...
```

Example for an `(a,e)` FLI map:

```text
# Arch of Chaos 1.6.0
a [AU]            e [-]             FLI
...
```

The output is written in scientific notation with high numerical precision.

---

# Automatic output-file naming

If the `-o` option is omitted, Arch of Chaos creates an informative name from the input parameters.

The numeric values preserve normal decimal points and scientific notation.

For example

```text
5.2026
1e-06
-30
```

remain recognizable in the generated file name.

---

## GRID naming

GRID names contain:

- indicator;
- formalism;
- every grid axis;
- lower limit of each grid axis;
- upper limit of each grid axis;
- number of intervals of each grid axis;
- integration-duration mode and value;
- `a2`;
- fixed orbital elements;
- fixed orbital phase, if the phase is not grid-controlled.

Example:

```text
LCI_NEWTONIAN_grid_a-4.8to5.8-N250_M-30to330-N30_nP-100_a2-5.2026_e-0_i-0_omega-60_Omega-0.txt
```

This corresponds to

```text
indicator = LCI
formalism = NEWTONIAN
a2 = 5.2026
nPeriods = 100

e = 0
i = 0
omega = 60
Omega = 0

grid = a 4.8 5.8 250
grid = M 30 330 30
```

---

## ORBIT / INDICATOR naming

ORBIT and INDICATOR use a common time-series naming scheme.

ORBIT example:

```text
ORBIT_NEWTONIAN_nP-100_dt-10_a2-5.2026_a-5.2_e-0.1_i-0_omega-60_Omega-0_M-30.txt
```

LCI example:

```text
LCI_NEWTONIAN_nP-100_dt-10_a2-5.2026_a-5.2_e-0.1_i-0_omega-60_Omega-0_M-30.txt
```

FLI example:

```text
FLI_HAMILTONIAN_T-1000_dt-10_a2-5.2026_a-5.2_e-0.1_i-0_omega-60_Omega-0_M-30.txt
```

If an explicit output is supplied,

```powershell
archofchaos.exe -i input.txt -o custom_name.txt
```

then `custom_name.txt` is used instead of an automatically generated name.

---

# Examples

## Example 1: ORBIT with `nPeriods`

```text
# ---------------------------------------------------------------------------
# Single-orbit integration
# ---------------------------------------------------------------------------

mode = ORBIT
indicator = NONE
formalism = NEWTONIAN

m1 = 1.0
m2 = 0.001001001001001
a2 = 5.2026

t0 = 0.0
nPeriods = 10

output_dt = 10.0

relTol = 1.0e-6
absTol = 1.0e-10

a = 5.2
e = 0.1
i = 0.0
omega = 60.0
Omega = 0.0
M = 30.0
```

Run:

```powershell
archofchaos.exe -i input_orbit.txt --verbose
```

---

## Example 2: LCI time series

```text
mode = INDICATOR
indicator = LCI
formalism = NEWTONIAN

m1 = 1.0
m2 = 0.001001001001001
a2 = 5.2026

t0 = 0.0
T = 1000.0

output_dt = 10.0

relTol = 1.0e-6
absTol = 1.0e-10

a = 5.2
e = 0.1
i = 0.0
omega = 60.0
Omega = 0.0
M = 30.0

dy1 = 1.0
dy2 = 0.0
dy3 = 0.0
dy4 = 0.0
```

---

## Example 3: FLI map in the `(a,e)` plane

```text
mode = GRID
indicator = FLI
formalism = NEWTONIAN

m1 = 1.0
m2 = 0.001001001001001
a2 = 5.2026

t0 = 0.0
nPeriods = 100

relTol = 1.0e-6
absTol = 1.0e-10

i = 0.0
omega = 0.0
Omega = 0.0
M = 60.0

grid = a 3.0 20.0 10000
grid = e 0.0 0.99 300

dy1 = 1.0
dy2 = 0.0
dy3 = 0.0
dy4 = 0.0
```

---

## Example 4: LCI map in the `(a,M)` plane

```text
mode = GRID
indicator = LCI
formalism = NEWTONIAN

m1 = 1.0
m2 = 0.001001001001001
a2 = 5.2026

t0 = 0.0
nPeriods = 100

relTol = 1.0e-6
absTol = 1.0e-10

e = 0.0
i = 0.0
omega = 60.0
Omega = 0.0

grid = a 4.8 5.8 250
grid = M 30.0 330.0 30

dy1 = 1.0
dy2 = 0.0
dy3 = 0.0
dy4 = 0.0
```

---

## Example 5: three-dimensional grid

The grid infrastructure is not restricted to two axes.

For example:

```text
grid = a 3.0 4.0 1
grid = e 0.0 0.2 2
grid = M 0.0 360.0 4
```

contains

\[
2\times3\times4=24
\]

grid points.

The first axis, `a`, varies fastest.

---

# Plotting GRID results

The repository can include the companion script

```text
plot_archofchaos_grid_heatmap.py
```

for plotting two-dimensional GRID outputs.

The plotting script:

- reads the original Arch of Chaos input file;
- reads the GRID result file;
- automatically identifies the two grid axes;
- reconstructs the rectangular grid;
- creates a heatmap;
- reads the indicator and formalism;
- reads `m1`, `m2`, `a2`;
- computes/displays the CRTBP mass parameter;
- reads the integration-duration mode;
- displays `T` or `nPeriods`;
- displays tolerances;
- displays fixed orbital elements;
- adds a compact information box to the figure;
- saves both PDF and PNG output;
- can use logarithmic indicator display;
- can automatically construct a descriptive plot file name.

Example:

```powershell
python .\plot_archofchaos_grid_heatmap.py `
    -i .\input_LCI_a_M.txt `
    -r .\LCI_a_M.txt
```

Optional explicit plot-output prefix:

```powershell
python .\plot_archofchaos_grid_heatmap.py `
    -i .\input_LCI_a_M.txt `
    -r .\LCI_a_M.txt `
    -o .\my_LCI_map
```

Logarithmic display:

```powershell
python .\plot_archofchaos_grid_heatmap.py `
    -i .\input_LCI_a_M.txt `
    -r .\LCI_a_M.txt `
    --log10
```

The plotting utility requires Python together with

```text
numpy
matplotlib
```

---

# Validation and regression tests

During development, several independent consistency checks have been performed.

These are currently development/regression checks rather than a formal automated unit-test framework.

## Integration-duration equivalence

ORBIT calculations were compared using equivalent durations specified as

```text
T
```

and

```text
nPeriods
```

and produced matching trajectories.

## LCI duration equivalence

LCI calculations using equivalent `T` and `nPeriods` input were compared and agreed.

## Newtonian versus Hamiltonian formulation

FLI calculations were compared between

```text
formalism = NEWTONIAN
```

and

```text
formalism = HAMILTONIAN
```

with matching results for the tested cases.

## Grid traversal tests

The generalized `GridIterator` has been tested with:

- two-dimensional `(a,e)` grids;
- reversed axis order `(e,a)`;
- three-dimensional `(a,e,M)` grids;
- periodic `M = 0 ... 360 deg`;
- partial `M = 0 ... 180 deg`;
- 10-degree mean-anomaly spacing;
- total-point counting;
- multidimensional carry propagation;
- first-axis-fastest ordering.

Example verified sequence:

```text
point 0 : a[0]=3  e[0]=0  M[0]=0
point 1 : a[1]=4  e[0]=0  M[0]=0
point 2 : a[0]=3  e[1]=0.1  M[0]=0
...
```

## Grid application tests

`GridIterator::apply()` has been tested to ensure that:

- grid-controlled `a` is applied;
- grid-controlled `e` is applied;
- angular values are converted from degrees to radians;
- non-grid-controlled elements remain unchanged.

## Generalized GRID integration

Small `(a,e)` FLI and LCI maps have been successfully computed with the generalized GRID implementation.

---

# Project structure

The exact repository layout may evolve, but the current code is organized around modules with responsibilities similar to the following:

```text
main.cpp
    command line
    run control
    output stream
    automatic output-file naming

init_data.h / init_data.cpp
    input parsing
    validation
    fixed orbital elements
    duration and phase input

grid.h / grid.cpp
    OrbitalElement
    GridAxis
    multidimensional GridIterator
    periodic endpoint handling
    grid application to orbital elements

model.h
    model/formalism/indicator abstractions

crtbp.h
    planar CRTBP dynamics and coordinate transformations

orbit.h / orbit.cpp
    orbital-element <-> Cartesian transformations

math_utils.h
    numerical/mathematical helper functions

version.h
    program name
    version information
    program metadata
```

The codebase uses Doxygen-style comments for public functions, classes, and important helpers.

---

# Known limitations

The current implementation has the following important limitations.

- The dynamical model is **planar**, not spatial.
- Inclination is therefore constrained to zero by the CRTBP2D validation.
- RLI is not yet implemented.
- GRID output is a flat table; multidimensional visualization is currently handled externally.
- The provided heatmap plotting utility is intended for exactly two grid axes.
- The project is currently developed primarily with Visual Studio/Windows.
- A fully automated unit/regression test suite is not yet part of the repository.
- Long high-resolution grids can require very large computation times.
- No parallel GRID execution is documented in the current version.

---

# Reproducibility recommendations

For research calculations, keep together:

```text
input file
output file
program version
compiler/build configuration
plotting script version
```

The automatically generated output file names are intended to make result directories more self-describing, but the **input file remains the authoritative record of the run configuration**.

For published results, it is recommended to record at least:

- Arch of Chaos version;
- mathematical formalism;
- masses and `a2`;
- grid definition;
- integration-duration convention;
- numerical tolerances;
- initial deviation vector;
- fixed orbital elements;
- compiler/build information.

---

# Suggested repository additions

For a public GitHub repository, the following files would be useful in addition to this README:

```text
LICENSE
CHANGELOG.md
CITATION.cff
CONTRIBUTING.md
examples/
scripts/
tests/
```

In particular:

- `LICENSE` should explicitly define redistribution/use conditions;
- `CHANGELOG.md` can summarize changes between versions;
- `CITATION.cff` can provide a standard citation entry for scientific use;
- `examples/` can contain small ORBIT, INDICATOR, and GRID inputs;
- `scripts/` can contain the Python heatmap utility;
- `tests/` can later host automated regression tests.

---

# Releases and version history

This `README.md` documents the current state of the repository rather than a single fixed release.

Release-specific changes are recorded in [`CHANGELOG.md`](CHANGELOG.md), while Git tags and GitHub Releases preserve the exact source and README corresponding to each published version.

Recommended release tags follow the form

```text
v1.5.0
v1.6.0
...
```

