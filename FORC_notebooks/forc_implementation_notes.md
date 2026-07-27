# FORC implementation notes

Design, mathematics, rationale, and caveats for the FORC processing in `pmagpy.forc`.

This page builds on the `FORCme_overview` and `FORCme_SOP` documents that accompanied the original FORCme release by Maxwell Brown (Institute for Rock Magnetism, University of Minnesota), work carried out under NSF award EAR-2148549. It describes the code as it now stands in RockmagPy, which differs from that release in several respects noted below.

## Purpose and scope

`pmagpy.forc` processes FORC data from raw instrument files through to diagrams and quantitative profiles. It is driven from a notebook through a single entry point, `process_forc`, which supports four input modes:

| `mode` | input | behaviour |
| --- | --- | --- |
| `"i"` | one file | process a single FORC measurement |
| `"s"` | a directory | stack repeated measurements of one specimen in magnetization space, then compute one distribution |
| `"b"` | a directory | process each matching file independently, returning one result per file |
| `"m"` | a MagIC measurements table | inspect specimen, experiment and sequence fields and dispatch automatically to `i`, `s`, or `b` |

Given one file or a directory, the pipeline:

1. resolves the file's units and parses the numeric section into calibration points and reversal curves;
2. corrects instrument drift using the interleaved calibration measurements;
3. optionally replaces the first and last point of each curve by linear extrapolation from interior points, since those can be affected by instrument settling;
4. optionally subtracts a reference reversal curve before the distribution calculation;
5. optionally regrids the curves onto a regular field grid;
6. builds a gridded *M*(H<sub>a</sub>, H<sub>b</sub>) array and estimates ρ(H<sub>a</sub>, H<sub>b</sub>) by local quadratic regression;
7. renders the result in both measurement coordinates and rotated B<sub>u</sub>–B<sub>c</sub> coordinates;
8. exports figures, profile data, and a MagIC measurements table.

## Field conventions

`Ha` is the **reversal field** at which a curve begins and `Hb` the **applied field** along that curve, with H<sub>b</sub> ≥ H<sub>a</sub>. This follows [Pike et al. (1999)](https://doi.org/10.1063/1.370176) and the subsequent literature.

> **Change from the original FORCme release.** That release used the opposite assignment, with `Hb` the reversal field and `Ha` the applied field, and inverted its rotation formulae to match. The computed distribution was correct, but every public array key meant the opposite of what the literature implies. The names were corrected before release in PmagPy; a consistent rename of every identifier is value-preserving, and ρ was verified bitwise unchanged across four datasets. Code written against the original release must swap `Ha` and `Hb`.

Fields are carried in tesla and moments in A m². MicroMag writes either SI or cgs depending on an instrument setting, and older Series 0015 exports default to cgs. Units are resolved from the file's `Units of measure` line or its column-unit line, header field tags are converted, and the converted field range is checked against what is physically plausible in a laboratory — a file that would have to reach tens of tesla is rejected rather than silently rescaled by 10<sup>4</sup>.

MicroMag header tags `Hb1`/`Hb2` and `Hc1`/`Hc2` are the **bias** and **coercivity** axis limits of the FORC *diagram*, unrelated to the reversal and applied fields. They are read into `Bu_max` and `Bc_max`.

## FORC mathematics

### Definition

Each FORC is measured by saturating the specimen, decreasing the field to a reversal field H<sub>a</sub>, then measuring the magnetization *M* as the applied field is increased from H<sub>a</sub> to the maximum field, with H<sub>b</sub> ≥ H<sub>a</sub>. The FORC distribution is

$$\rho(H_a, H_b) = -\frac{1}{2}\,\frac{\partial^2 M(H_a, H_b)}{\partial H_a \, \partial H_b}.$$

In practice *M* is known only on a discrete, and often irregular, set of points, and must be smoothed before numerical differentiation.

### Rotation to (B<sub>u</sub>, B<sub>c</sub>)

$$B_u = \frac{H_b + H_a}{2}, \qquad B_c = \frac{H_b - H_a}{2}.$$

### Local quadratic estimate

At each grid node a weighted quadratic surface is fitted to the neighbouring measurements,

$$M \approx \beta_0 + \beta_1 \Delta H_b + \beta_2 \Delta H_a + \beta_3 (\Delta H_b)^2 + \beta_4 (\Delta H_b)(\Delta H_a) + \beta_5 (\Delta H_a)^2,$$

and the mixed term gives the distribution,

$$\rho(H_a, H_b) = -\tfrac{1}{2}\beta_4.$$

Neighbours are selected within an ellipse in index space. Each neighbour is assigned a distance proxy *u* ∈ [0, 1] and the standard LOESS tricube weight ([Cleveland, 1979](https://doi.org/10.1080/01621459.1979.10481038)),

$$w(u) = (1 - u^3)^3.$$

Local coordinates are formed from the difference of the measured field values, scaled by the nominal field step before the normal equations are assembled. The scaling keeps the 6 × 6 systems well conditioned; taking the differences from the measured values rather than from index offsets is what makes the estimate correct on the unevenly spaced grids that real measurements produce. The mixed coefficient is converted back to inverse tesla squared afterwards. Cells whose local system is rank deficient or near singular are rejected before the solve, so a node with too few or collinear neighbours returns NaN rather than a meaningless coefficient.

A quadratic is the minimal model that captures the cross-term needed for the mixed derivative, and is the common choice across FORC smoothing implementations.

## Adaptive smoothing

The LOESS neighbourhood is controlled by spans in H<sub>a</sub> and H<sub>b</sub> (`span_Ha_T`, `span_Hb_T`) and a minimum number of points (`min_pts`). To reduce guesswork, `guess_loess_params` estimates these from the geometry of the measurement:

1. estimate the grid steps *dH*<sub>a</sub> and *dH*<sub>b</sub> from medians of differences;
2. compute the fill fraction *f* = *N*<sub>finite</sub> / *N*<sub>physical</sub> of the gridded *M* array;
3. for a candidate ellipse (r<sub>x</sub>, r<sub>y</sub>) in index units, count the offsets inside the ellipse (*N*<sub>cand</sub>) and estimate the expected number of usable points as *N*<sub>eff</sub> ≈ *f* *N*<sub>cand</sub>;
4. grow the ellipse until *N*<sub>eff</sub> meets a target (`target_n_eff`, default 60), subject to caps;
5. set `min_pts` to a fraction of *N*<sub>eff</sub> (`min_pts_frac`, default 0.55).

The user then scales these through `smooth_strength` and `min_pts_strength`, and every value actually used is recorded in `out['loess_params']`.

> **Change from the original FORCme release.** The fill fraction was computed over the whole rectangular array, roughly half of which is the structurally empty region H<sub>b</sub> < H<sub>a</sub>. That understated the local data density and inflated the smoothing window. It is now computed over the physical half-plane only. For the Baraboo example this changes the automated span from 35 mT to 25 mT; the peak position is unchanged and the peak amplitude rises by about 12%.

This approach keeps the local regression well determined without oversmoothing, but it operates at the level of the dataset rather than varying across the diagram. It is a starting point, and results should be assessed across a range of smoothing levels, as one would do manually in other FORC software. Where a single window cannot serve the whole diagram, the spatially variable scheme of Egli (2013) is implemented and is selected with `smoothing='variforc'`.

## Drift correction

A calibration measurement is made at a fixed field before each reversal curve. The change in the measured moment at that field over the run is instrument drift, and it is removed by interpolating the calibration record across the measurement sequence, either linearly or with a monotone piecewise-cubic scheme (`drift_fit`).

Two caveats: the drift abscissa is the segment index rather than elapsed time, which slightly misweights curves of unequal length; and drift accrued *within* a single sweep is not modelled.

> **Change from the original FORCme release.** Calibration points that occupy their own one-row block — the layout MicroMag produces when every reading is separated by a blank line, which includes the Baraboo example — were discarded before they could be classified, so no drift correction was applied at all and nothing said so. Isolated points at the calibration field are now retained and classified as calibration measurements, and the pipeline reports explicitly when fewer than two are found rather than silently returning uncorrected curves. `out['drift_corrected']` and `out['n_calibration_points']` record what happened.

## Reference-curve subtraction

Subtracting a reference reversal curve from the whole family removes the common reversible and high-field signal and can aid visual inspection of the curves. It does not change the distribution: the subtracted baseline is a function of the applied field alone and vanishes under the mixed derivative.

`reference_curve="lowest_reversal"` (the default) selects the curve measured from the most negative reversal field, which spans the widest field range and is the closest measured approximation to the lower branch of the major loop. `"first_measured"` selects the first curve in measurement order.

> **Change from the original FORCme release.** The option was named `do_lower_branch_subtract` and always used the first-measured curve. For MicroMag files that is the curve at the *highest* reversal field, typically only a few points long and spanning about 10 mT, so the interpolated "lower branch" was clamped to a near-constant offset over most of the field range. The option is now `do_reference_subtract` with an explicit `reference_curve` choice, and interpolation outside the reference curve's own span yields NaN rather than a clamped endpoint value.

## Gridding

`build_forc_grid` assembles the measured curves into a rectangular *M*(H<sub>a</sub>, H<sub>b</sub>) array. Where the measured applied fields already lie on the common lattice, each measurement is placed in its own column exactly. Where they do not — which happens whenever the reversal-field increment is not a whole multiple of the applied-field step, so that every curve starts at a different sub-step offset — the curves are interpolated onto the lattice instead.

> **Change from the original FORCme release.** Measurements were always assigned to the nearest column. For an incommensurate file that displaces measurements by up to half a field step, and because the distribution is a second derivative in field, the resulting error is not small: on a synthetic surface with a known answer it changed ρ by a factor of about three. Interpolating instead recovers the analytic value exactly. Setting `do_regrid=True` also avoids the problem and remains a good default for consistency across specimens.

## Display and plotting

The B<sub>u</sub>–B<sub>c</sub> plot uses rotated cell **edges** for the colour image and rotated cell **centres** for contours. Drawing each H<sub>a</sub>–H<sub>b</sub> cell as a true parallelogram avoids the coercivity-axis offset that arises when a rotated centre grid is passed directly to `pcolormesh`, and lets the image reach the B<sub>c</sub> = 0 axis correctly.

`display_upsample_factor` applies a plotting-only refinement of the grid; the default is 0, the rawest representation of the data, which can look blocky if a coarse field increment was used. Upsampling increases in steps of 2, and 2 is usually sufficient. It trades against data area, most noticeably near B<sub>u</sub> = 0. Colour normalization is computed on the native grid before any upsampling, so upsampling cannot alter the reported normalization.

`normalize_to_unit=True` scales by a percentile of |ρ| so the colour scale runs from −1 to +1; values above that percentile clip. This affects only the visualization, not the underlying ρ values, which are always returned in their calculated units.

## Profiles

One-dimensional cuts through ρ quantify the position, amplitude, and width of features more precisely than can be read from the colour plot alone. An additional one-dimensional smoothing, `smooth_sigma_bins`, may be applied to the extracted profiles; it affects only the profile curves and does not recalculate the two-dimensional distribution.

- **Automated profiles** (`plot_auto_forc_profiles`) locate the bounded peak of ρ within the displayed window and construct a horizontal profile, ρ against B<sub>c</sub> at the peak B<sub>u</sub>, and a vertical profile, ρ against B<sub>u</sub> at the peak B<sub>c</sub>.
- **Custom profiles** (`plot_custom_forc_profiles`) take user-specified B<sub>u</sub> and B<sub>c</sub> values, for interrogating a shoulder, a secondary ridge, or a low-coercivity contribution that is scientifically important but not dominant in amplitude.
- **Ridge-following profiles** (`ridge_profile`) track the B<sub>u</sub> of the maximum of ρ in each coercivity bin and evaluate ρ along the crest. Where the crest drifts in B<sub>u</sub> with coercivity — hematite is a common case — a fixed-B<sub>u</sub> cut slides off the ridge and understates the coercivity distribution.

Both the numerical profile data and the profile figures are exported, so the record shows exactly which cuts were taken and how they were displayed.

> **Changes from the original FORCme release.** Profile smoothing in the custom-profile path called a function that did not exist; the resulting `NameError` was caught by a bare `except` and the profile returned unsmoothed regardless of the requested `smooth_sigma_bins`. Profiles were also binned as a fixed count over the full measured field range rather than over the requested window, which quantized the reported peak position and diluted the peak amplitude for a measurement whose reversal curves extend far below the plotted region. Bins now follow the requested window at the grid step by default, settable through `bin_width`. Note that `smooth_sigma_bins` is a Gaussian sigma **in bins**, so its physical width scales with the bin width.

## Output locations

Outputs are written beneath the input directory:

- `FORC_figures` — main FORC figures.
- `FORC_profiles` — profile text files and profile figures.
- `MagIC` — the exported MagIC measurements table, when writing from a raw instrument file.

Titles default to the file name without its extension in modes `i` and `b`, and to the common part of the matched file names in mode `s`. `sample_title` overrides this.

## Recommended settings

- Use `stack_method="median"` when repeated files may contain outliers.
- Start with `display_upsample_factor=0`; increase to 2 or 4 only for smoother plotting.
- Use `do_regrid=True` for consistency across specimens, or if the source export is not already on a regular field grid.
- Use reference-curve subtraction only when justified by the measurement protocol, and remember it does not change ρ.
- Assess the result across a range of `smooth_strength` values before settling on one, and report the value used.

## Troubleshooting

- **No files matched** — check `path` and `file_type`.
- **Too many NaNs in ρ** — reduce `min_pts_strength`, or inspect data sparsity and the reported segmentation.
- **"NO drift correction has been applied"** — fewer than two calibration points were identified. Check that the file contains calibration measurements and that `HCal` in the header matches the field at which they were made; `cal_tol_T` sets the matching tolerance.
- **"not a plausible laboratory field"** — the file is almost certainly in cgs but does not declare it. Check the `Units of measure` line.
- **Unexpected stack title** — set `sample_title` explicitly.
- **Plots too blocky or too smooth** — adjust `display_upsample_factor`, which is a display setting, before adjusting `smooth_strength`, which is not.

## Remaining caveats

- The neighbourhood ellipse and its tricube weights are defined in index space, while the polynomial basis uses true field distances. These are consistent only for uniformly spaced grids; for adaptive or even-moment FORC files with non-uniform reversal-field spacing, regridding is recommended.
- The LOESS neighbourhood is truncated at the H<sub>b</sub> = H<sub>a</sub> diagonal, so values near B<sub>c</sub> = 0 are estimated from a half-window. `edge_mask_bc_bins` can mask that region; reflection in the manner of VARIFORC is not implemented.
- The default LOESS estimator adapts its window at the dataset level but not spatially, and propagates no uncertainty. Both are addressed by the VARIFORC option documented below, which is opt-in.

---

---

# VARIFORC variable smoothing

`pmagpy.forc` implements the variable-smoothing protocol of [Egli (2013)](https://doi.org/10.1016/j.gloplacha.2013.08.003), together with the per-cell error analysis of [Heslop and Roberts (2012)](https://doi.org/10.1029/2012GC004115). It is selected with `smoothing='variforc'` and is off by default, so nothing changes for existing analyses.

The implementation is independent: it was written from the equations published in Egli (2013) and the parameter semantics documented in the VARIFORC manual, not from the VARIFORC source, which is separately licensed. Where the two describe the same quantity differently, the paper was followed.

## Why a single smoothing window is not enough

Conventional processing applies one regression window across the whole diagram. That is a compromise whenever a diagram carries features on different scales at once, which is the normal case: a narrow ridge near B<sub>u</sub> = 0 sitting on a broad, low-amplitude background.

The compromise is quantifiable. On a synthetic ridge-plus-background surface whose distribution is known exactly, the window that minimizes error over the ridge and the window that minimizes noise over the background differ by roughly a factor of five, and no single choice serves both. [FORC_smoothing.ipynb](./FORC_smoothing.ipynb) works this through.

Egli's observation is that Preisach-type FORC functions become smoother with distance from the origin of the diagram, so the window can be allowed to grow with B<sub>c</sub> and |B<sub>u</sub>| without discarding anything real, while being held narrow across the features that need resolving.

## What the method does

**Regression in rotated coordinates.** The local quadratic is fitted over an upright rectangle in (B<sub>c</sub>, B<sub>u</sub>), not in the measurement coordinates. Only in the rotated frame can the smoothing be made anisotropic in the way the method needs — fine across a ridge, coarse along it. A rectangle upright in (B<sub>c</sub>, B<sub>u</sub>) is a 45° diamond in (H<sub>a</sub>, H<sub>b</sub>).

**A different derivative.** In rotated coordinates the distribution is not a mixed derivative. Egli's Eq. 1 gives

$$\rho = -\frac{1}{2}\frac{\partial^2 M}{\partial H_a \partial H_b} = \frac{1}{8}\left(\frac{\partial^2 M}{\partial B_c^2} - \frac{\partial^2 M}{\partial B_u^2}\right),$$

so ρ is read from the *difference of the two pure second-order coefficients* of the local fit rather than from the cross term.

**A weight with a unit core.** The regression weight is separable, `T(ΔBc/δH, s_c)·T(ΔBu/δH, s_b)`, with (Egli Eq. 6)

$$T(u,s) = \begin{cases} 1 & |u| \le s-1 \\ 1 - 2(|u|-s+1)^2 & s-1 < |u| \le s-\tfrac12 \\ 2(|u|-s)^2 & s-\tfrac12 < |u| \le s \\ 0 & |u| > s \end{cases}$$

Compared with the tricube weight of FORCinel this keeps peripheral measurements at full weight for longer, so fewer points are wasted for a given resolution, while the two quadratic pieces still join smoothly enough that the distribution is continuous as the window moves.

**Factors that grow, except across ridges.** The smoothing factors follow Egli Eq. 11–12 and 14,

$$s = \min\left[(1-\lambda)\,s_1 + \lambda\frac{|H|}{\delta H},\ \max\left(s_0, \frac{|H|}{\delta H}\right)\right].$$

The first term is the linear growth; the second holds the factor at a floor `s0` across a ridge and then relaxes just fast enough that a window centred near the ridge never reaches across it. Setting λ = 0 and `s1 = SF + 1` recovers conventional processing at the smoothing factor of Pike et al. (1999).

**Limits.** Growing windows generate their own artefacts where the magnetization surface changes sharply. Ridges away from the axes, and the two diagonals through the distribution maximum, can be given explicit caps; the diagonal caps can be placed automatically from a preliminary pass.

## Using it

The parameters of the original protocol are expressed in smoothing factors and growth rates. `variforc_settings` lets a diagram be described instead by what is visible in it, and converts that to the factors the kernel needs.

```python
import pmagpy.forc as forc

out = forc.process_forc(
    mode='m', path='measurements.txt',
    smoothing='variforc',
    variforc=forc.variforc_settings(
        'central_ridge',            # what the sample shows
        smoothing_factor=9,         # baseline window away from the ridge
        growth_rate=0.1,            # how fast it grows with distance
        central_ridge=4,            # window held this narrow across the ridge
        central_ridge_position=0.0004,   # the ridge sits 0.4 mT above Bu = 0
    ),
)
```

The four presets correspond to the cases Egli distinguishes:

| preset | use when the diagram shows |
| --- | --- |
| `regular` | no sharp ridges; a smooth Preisach-type distribution |
| `central_ridge` | a horizontal ridge along B<sub>u</sub> = 0 (non-interacting single-domain particles, magnetofossils) |
| `vertical_ridge` | a vertical ridge along B<sub>c</sub> = 0 (thermal relaxation, superparamagnetism) |
| `both_ridges` | both, over a continuous background (many volcanic ash and pseudo-single-domain samples) |

Everything used is recorded in `out['smoothing_params']`, so a published diagram can be reported with the parameters that produced it.

`growth_rate` between about 0.06 and 0.12 is sensible; 0 gives conventional processing, and above roughly 0.2 the windows become large enough to need the diagonal limits.

## Uncertainty and significance

With `estimate_uncertainty=True` the measurement noise is propagated through each local fit to a standard error on ρ. For a weighted least-squares fit the covariance of the coefficients is σ²N⁻¹SN⁻¹ with N = XᵀWX and S = XᵀW²X, and ρ is a fixed linear functional of those coefficients, so the standard error costs one extra contraction per batch rather than a second pass over the data.

The noise level itself is estimated from fourth-order differences along each reversal curve, which annihilate any local cubic and so respond to noise rather than to real curvature, with a median-absolute-deviation estimator so that sharp structure in part of the diagram does not inflate it.

The result adds `rho_sigma` and `snr` to the output. Following Egli's practice, SNR > 3 marks the part of the diagram that differs significantly from zero, and contours are only worth drawing where the SNR is well above that.

## How it is tested

The implementation is checked on three independent levels.

**From first principles.** On an analytic quadratic surface the rotated second-derivative-difference formulation and the mixed-derivative LOESS estimator must both be exact; they agree with each other and with the analytic value to 5 × 10⁻¹¹. On a separable Preisach assemblage with a closed-form distribution the peak is recovered at its known coordinates and the amplitude to within 5% across the interior.

**Against the documented algorithm.** Egli's Eq. 7 defines a synthetic unit ridge, `M = -4|B_u|`, whose processed resolution his Eq. 9 predicts as ΔH<sub>b</sub> = (1.076 s<sub>b</sub> − 0.468) δH. The implementation reproduces this without having been fitted to it:

| s<sub>b</sub> | measured FWHM/δH | Egli Eq. 9 |
| --- | --- | --- |
| 2 | 1.678 | 1.684 |
| 3 | 2.746 | 2.760 |
| 5 | 4.898 | 4.912 |
| 10 | 10.295 | 10.292 |

A least-squares fit to the measured values gives 1.078 s<sub>b</sub> − 0.486 against the published 1.076 s<sub>b</sub> − 0.468.

**Through invariances.** Nodes are grouped internally only to share a candidate offset list; each is then weighted by its own exact smoothing factors. The result is therefore independent of the grouping, and the test suite asserts that it is, to 10⁻⁶ of the peak. This matters: an earlier design quantized the smoothing factors themselves, which was faster but printed visible bands across the diagram at a few percent of peak amplitude.

**Uncertainty.** The propagated standard error is checked against the actual scatter of ρ over independent noise realizations; the median ratio of predicted to observed is 1.02.

**On the published examples.** Every VariFORC example dataset that ships with measurements is processed with Egli's own parameter files. See [FORC_smoothing.ipynb](./FORC_smoothing.ipynb) for the diagrams.

**Against an independent implementation.** The variable-smoothing kernel was cross-checked cell by cell against FORCsensei ([Heslop et al., 2020](https://doi.org/10.1029/2020JB020418), MIT licensed), which implements the same Egli protocol from an independent codebase. Both of its case-study datasets were processed through both packages under constant and variable smoothing:

| dataset | grid | smoothing | max \|Δρ\|/ρ<sub>max</sub> | median |
| --- | --- | --- | --- | --- |
| UB199BF | 136 × 108 | s = 4 | 1.7 × 10⁻¹² | 1.2 × 10⁻¹⁴ |
| UB199BF | 136 × 108 | s<sub>0</sub> = 3, s<sub>1</sub> = 7, λ = 0.1 | 3.8 × 10⁻¹² | 1.2 × 10⁻¹⁴ |
| MD00-2361 | 292 × 175 | s = 4 | 4.4 × 10⁻¹² | 2.6 × 10⁻¹⁴ |
| MD00-2361 | 292 × 175 | s<sub>0</sub> = 3, s<sub>1</sub> = 7, λ = 0.1 | 1.1 × 10⁻¹¹ | 2.0 × 10⁻¹⁴ |

Agreement is at the level of floating-point accumulation order, for both the NumPy and the numba engine. Reaching it required fixing a real defect that none of the internal tests could have caught, because every synthetic fixture used a perfectly uniform grid: the kernel had been forming each neighbour's local coordinates as *index offset × median field step*. Measured reversal fields are not exactly evenly spaced, so on real data both the smoothing weights and the polynomial basis were evaluated at slightly displaced positions. Median disagreement with FORCsensei was 8 × 10⁻⁵ on real data against 10⁻¹³ on synthetic data — a discrepancy that pointed directly at an assumption the synthetic tests shared with the code. Coordinates are now taken from the true field values, which makes the design matrix per-node rather than shared.

The same uniform-spacing assumption appeared a second time, in the pre-filter that trimmed the candidate offset list to those falling inside a node's rectangle. Judged with median spacing, that test can exclude an offset that the true spacing places inside the window, silently truncating that node's neighbourhood. The pre-filter was removed; points outside a node's own rectangle receive zero weight, which is exact and costs roughly a factor of two in candidates.

## Planned: central-ridge isolation

Not yet implemented. This is the most substantial remaining gap against FORCinel
(`IsolateCR`) and VARIFORC, and it is the feature that makes variable smoothing worth
having in the first place, so it should be built next. Notes toward an implementation.

### What the central ridge is, and why it needs separating

A non-interacting uniaxial single-domain assemblage switches at a well-defined field with
no bias, so its FORC signature is a ridge of zero width in B<sub>u</sub>, running along
B<sub>u</sub> = 0 ([Egli et al., 2010](https://doi.org/10.1029/2009GC002916)). Any
measured ridge therefore has the width the *processing* gave it, not a width the sample
has. Magnetostatic interactions, thermal relaxation, and vortex or multidomain particles
all produce a genuinely two-dimensional background that spreads in B<sub>u</sub>.

Separating the two matters because the ridge integral is proportional to the abundance of
the non-interacting single-domain fraction, which is the basis for quantifying
magnetofossil content in sediments ([Egli et al., 2010](https://doi.org/10.1029/2009GC002916);
[Heslop et al., 2014](https://doi.org/10.1002/2014GC005291)). Reading an amplitude off the
unseparated diagram conflates the ridge with whatever background sits under it.

### The approach to implement

Work column by column in B<sub>c</sub>. For each coercivity bin, the profile
&rho;(B<sub>u</sub>) at fixed B<sub>c</sub> is modelled as a narrow peak centred on
B<sub>u</sub> = 0 sitting on a smooth background:

1. **Estimate the background** from the parts of the profile outside the ridge. Fit a low
   order polynomial in B<sub>u</sub> to the flanks, excluding a window around
   B<sub>u</sub> = 0 whose width is set by the smoothing actually applied, not by a fixed
   number. Egli's Eq. 9 already gives that width, &Delta;H<sub>b</sub> =
   (1.076 s<sub>b</sub> &minus; 0.468) &delta;H, and `variforc_smoothing_factors` returns
   the per-node s<sub>b</sub>, so the exclusion window can be derived rather than guessed.
   This is the single most important design point: the ridge is *defined* as the component
   narrower than the processing resolution, so the resolution must drive the separation.
2. **Subtract**, leaving a residual that should be a peak near B<sub>u</sub> = 0.
3. **Integrate the residual over B<sub>u</sub>** to give the ridge profile
   f<sub>CR</sub>(B<sub>c</sub>), and keep the background as its own two-dimensional array
   so both components can be plotted and both integrate to interpretable quantities.

### Cautions

- **Do not fit a fixed-width Gaussian to the ridge.** The measured ridge width is set by
  the smoothing, and with VARIFORC that width varies across the diagram. A fixed width
  will over- or under-subtract systematically as a function of B<sub>c</sub>.
- **The ridge is not always at B<sub>u</sub> = 0.** Weak mean interactions displace it, and
  the Baraboo hematite shows a crest that drifts to &minus;98 mT with increasing
  B<sub>c</sub>. `track_bu_offset_vs_bc` already follows the crest and should feed the
  centring rather than assuming zero. Egli's mean-field correction addresses the same
  problem from the other end and is also unimplemented.
- **Near B<sub>c</sub> = 0 the separation is ill-posed**, because the ridge and the
  background both pile up against the boundary and the usable half-window shrinks. Return
  NaN rather than a number there, as the smoothing code already does.
- **Propagate the uncertainty.** The per-cell standard error is already available; the
  ridge integral should carry one so that a reported abundance has an error bar. A ridge
  integral without one is exactly the kind of number that gets over-interpreted.

### Validation plan

The VariFORC example datasets include the magnetofossil and magnetotactic bacteria samples
that Egli used to demonstrate the method, and they are already redistributed under
`example_data/FORC/variforc_examples/`. Check in this order:

1. A synthetic ridge of known integrated amplitude on a known background: the separation
   must recover both to within the propagated uncertainty.
2. The magnetofossil sample, where the ridge is strong and its published behaviour is
   documented.
3. Invariance checks: the recovered ridge integral must be insensitive to the smoothing
   level over a reasonable range. If it is not, the exclusion window is wrong.

Only step 1 is a real test of correctness; steps 2 and 3 test that it behaves sensibly on
data. Do not treat visual similarity to a published figure as validation.

## Cost

The calculation is vectorized over grid nodes, in batches sharing a candidate offset list. For Egli's magnetofossil example — 449 reversal curves on a 0.54 mT grid, 160,000 output nodes — conventional processing takes about 1 s and variable smoothing 20 to 35 s depending on the growth rate. That is slower than the constant-window estimator, as it must be, since the windows are several times larger, but fast enough to iterate on.

## Differences from VARIFORC

- The output is computed on the measurement grid, optionally regridded, rather than on an independent output mesh. VARIFORC's mesh is decoupled from the measurement protocol, which lets datasets taken with different protocols be combined without interpolation.
- Egli's INPUT 09 offers three weight-margin widths; the weight here is the single form given as Eq. 6 in the paper, corresponding to the recommended middle setting.
- The lower-diagonal trim is implemented as a simple distance criterion rather than VARIFORC's trim factor.
- Central-ridge extraction (`IsolateCR`) and mean-field correction are not implemented.

## References

- Cleveland, W. S. (1979). Robust locally weighted regression and smoothing scatterplots. *Journal of the American Statistical Association*, 74, 829–836. https://doi.org/10.1080/01621459.1979.10481038
- Egli, R. (2013). VARIFORC: An optimized protocol for calculating non-regular first-order reversal curve (FORC) diagrams. *Global and Planetary Change*, 110, 302–320. https://doi.org/10.1016/j.gloplacha.2013.08.003
- Harrison, R. J., & Feinberg, J. M. (2008). FORCinel: An improved algorithm for calculating first-order reversal curve distributions using locally weighted regression smoothing. *Geochemistry, Geophysics, Geosystems*, 9, Q05016. https://doi.org/10.1029/2008GC001987
- Pike, C. R., Roberts, A. P., & Verosub, K. L. (1999). Characterizing interactions in fine magnetic particle systems using first order reversal curves. *Journal of Applied Physics*, 85, 6660–6667. https://doi.org/10.1063/1.370176
- Roberts, A. P., Pike, C. R., & Verosub, K. L. (2000). First-order reversal curve diagrams: A new tool for characterizing the magnetic properties of natural samples. *Journal of Geophysical Research*, 105, 28461–28475. https://doi.org/10.1029/2000JB900326
- Roberts, A. P., Heslop, D., Zhao, X., & Pike, C. R. (2014). Understanding fine magnetic particle systems through use of first-order reversal curve diagrams. *Reviews of Geophysics*, 52, 557–602. https://doi.org/10.1002/2014RG000462
