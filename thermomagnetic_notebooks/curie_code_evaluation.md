# Evaluation of Curie temperature code in PmagPy

This page evaluates the Curie temperature estimation code that exists across the PmagPy
codebase, documents what the RockmagPy notebooks use going forward, and makes
recommendations on what to retain or sunset. It reflects the state of the codebase as of
July 2026.

## What to use

The functions in `pmagpy/rockmag.py` are the supported path for Curie temperature
estimation:

- **`rockmag.curie_temperature_estimates`** — the main entry point. Takes a
  MagIC-formatted experiment DataFrame, applies the requested methods to the heating and
  cooling branches, and returns a tidy comparison table with per-method caveats attached.
- **Per-method functions on plain arrays** — `curie_derivative_estimates` (inflection
  point and maximum curvature), `curie_two_tangent`, `curie_inverse_susceptibility`
  (Curie-Weiss extrapolation with fit statistics and uncertainty), and
  `curie_landau_fit` (the in-field Landau equation-of-state fit of
  [Fabian et al., 2013](https://doi.org/10.1029/2012GC004440), with formal uncertainty
  and a Moskowitz-style extrapolation mode). These operate on temperature/magnetization
  arrays from any source, so instrument exports that are not in MagIC format can be
  analyzed directly.
- **`rockmag.plot_curie_estimates`** (static comparison figure),
  **`rockmag.interactive_curie_inverse_susceptibility`** (draggable Curie-Weiss fit for
  exploration), and **`rockmag.add_curie_estimates_to_specimens_table`** (archives the
  estimate, method, and uncertainty in a MagIC specimens table via `critical_temp` /
  `critical_temp_type`).

Method selection guidance — which estimators are appropriate for strong-field *Mₛ(T)*
versus low-field *χ-T* data, and the systematic biases of each — is covered in the
[Curie temperature estimation notebook](./curie_temperature_estimation.ipynb) and in
the function docstrings. These functions are covered by a test suite
(`pmagpy/test/test_rockmag_curie.py`) with synthetic curves that have known analytical
Curie temperatures, including tests that the documented bias ordering (two-tangent ≈
maximum curvature > inflection point on in-field curves) is reproduced.

## Legacy inventory and critique

### `ipmag.curie` (pmagpy/ipmag.py)

The classic PmagPy implementation (contributed 2008): reads a two-column T,M file (or a
MagIC measurements file), resamples the curve onto a uniform 1 °C grid, smooths with a
Bartlett-window convolution (`ipmag.smooth`), forms first and second derivatives by
smoothed central differences, and reports the Curie temperature as the maximum of the
second derivative, together with a sliding-window stability plot.

Limitations, in light of the current literature:

1. **Single method, with a known bias.** The second-derivative maximum is the
   maximum-curvature estimate, which on in-field *M(T)* curves lies systematically above
   the Curie temperature — typically by 10–15 °C, increasing with applied field
   ([Fabian et al., 2013](https://doi.org/10.1029/2012GC004440)). No alternative
   estimators are offered, and no uncertainty is reported.
2. **Fixed-grid resampling.** The interpolation onto a 1 °C grid requires integer-degree
   temperature steps (the function exits otherwise) and resamples sharp features.
3. **No heating/cooling handling.** The function treats the input as a single curve;
   heating-cooling irreversibility — the primary in-run diagnostic of alteration — is
   not addressed.
4. **Output coupled to plotting.** The estimate is printed and plotted rather than
   returned in a form suitable for tabulation across specimens.

For continuity: on the classic `curie_example.dat` file, `ipmag.curie` with a 10-degree
window yields 552 °C; the `max_curvature` estimate from
`rockmag.curie_derivative_estimates` reproduces this value within a few degrees (the
small difference reflects the different smoothing kernels), so legacy results map
directly onto the `max_curvature` column of the new comparison table.

### `programs/curie.py` (command-line program)

The CLI duplicates the `ipmag.curie` algorithm with its own local copy of the smoothing
function rather than importing it — two implementations of the same algorithm to
maintain.

### `ipmag.smooth`

A uniform-grid convolution smoother used by `ipmag.curie`. Superseded for thermomagnetic
work by `rockmag.smooth_moving_avg`, which operates in temperature space (no resampling
required, tolerant of irregular spacing), supports the same window families, and can
return within-window variances used by `rockmag.optimize_moving_average_window` for
objective window selection.

### `rockmag.estimate_curie_temperature` (removed)

An earlier rockmag function (added in PmagPy PR #776) that returned a six-tuple of
first-derivative-minimum, second-derivative-maximum, and zero-crossing temperatures for
the heating and cooling branches. It has been replaced by
`curie_temperature_estimates`, which returns the same derivative-based estimates (as
the `inflection` and `max_curvature` methods) in a tidy DataFrame, adds the
inverse-susceptibility and Landau-fit methods, and attaches the method caveats to the
output. No notebooks or other code called the old function.

## Relationship to PmagPy PR #856

The open PR [#856](https://github.com/PmagPy/PmagPy/pull/856) ("Removal of rock mag
functions from ipmag") deletes `ipmag.curie` and rewrites `programs/curie.py` as a thin
wrapper around `rockmag.plot_ms_t`. That direction — consolidating rock magnetic
analysis in `rockmag.py` — is the right one, and this evaluation supports it. However,
as currently written the PR leaves the CLI with plotting only: the quantitative Curie
temperature estimate is lost.

**Recommendation:** once the `rockmag` Curie functions land, `programs/curie.py` should
wrap `rockmag.curie_temperature_estimates`, with the `max_curvature` method as the
default (reproducing legacy behavior for existing users) and a `--method` flag exposing
`inflection`, `two_tangent`, `inverse_susceptibility`, and `landau`. This restores the
CLI's quantitative output while retiring the duplicated legacy algorithm. No changes to
the PR's branch are made as part of this evaluation.

## Summary of recommendations

| code | status | recommendation |
|---|---|---|
| `rockmag.curie_temperature_estimates` + per-method functions | current | use for all new work |
| `rockmag.estimate_curie_temperature` | removed | superseded by the above |
| `ipmag.curie` | legacy | sunset via PR #856; results map to the `max_curvature` method |
| `programs/curie.py` | legacy | rewrap around `rockmag.curie_temperature_estimates` with `--method` (recommended amendment to PR #856) |
| `ipmag.smooth` | legacy | superseded by `rockmag.smooth_moving_avg`; can be removed with `ipmag.curie` |
