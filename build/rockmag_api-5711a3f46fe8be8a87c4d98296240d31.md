# `pmagpy.rockmag` API Reference

## `ANOVA`

```python
ANOVA(xs, ys)
```

ANOVA statistics for linear regression

**Parameters**

- **xs** (`numpy array`) — x values
- **ys** (`numpy array`) — y values

**Returns**

- `dict` — dictionary of the results of the ANOVA calculation and intermediate statistics for the ANOVA calculation

---

## `Fabian_nonlinear_fit`

```python
Fabian_nonlinear_fit(H, chi_HF, Ms, alpha, beta)
```

function for calculating the Fabian non-linear fit

**Parameters**

- **H** (`numpy array`) — field values
- **chi_HF** (`float`) — high field susceptibility
- **Ms** (`float`) — saturation magnetization
- **alpha** (`float`) — coefficient of the H^(beta) term, needs to be negative
- **beta** (`float`) — exponent of the H^(beta) term, needs to be negative

**Returns**

- `numpy array` — fitted magnetization values for each field value in `H` (same shape as `H`)

---

## `IRM_nonlinear_fit`

```python
IRM_nonlinear_fit(H, chi_HF, Ms, a_1, a_2)
```

Calculate the non-linear fit for Isothermal Remanent Magnetization (IRM) as a function of applied field.

This function models the IRM signal as a sum of high-field linear susceptibility, 
saturation magnetization, and non-linear correction terms with inverse field dependence.
The model is commonly used for fitting high-field IRM data, especially for extracting 
parameters such as high-field susceptibility (chi_HF) and saturation magnetization (Ms).

**Parameters**

- **H** (`numpy.ndarray`) — Array of applied magnetic field values (in Tesla).
- **chi_HF** (`float`) — High-field magnetic susceptibility. Converted to Tesla to match the unit of the field.
- **Ms** (`float`) — Saturation magnetization (in the same units as IRM).
- **a_1** (`float`) — Coefficient for the H^(-1) non-linear correction term. Should be negative.
- **a_2** (`float`) — Coefficient for the H^(-2) non-linear correction term. Should be negative.

**Returns**

- `numpy.ndarray` — Array of fitted IRM values corresponding to each field value in `H`. 

**Examples**

```python
>>> H = np.linspace(0.1, 3, 100)  # field in Tesla, avoid zero for stability
>>> fit = IRM_nonlinear_fit(H, chi_HF=0.02, Ms=1.2, a_1=-0.03, a_2=-0.01)
>>> import matplotlib.pyplot as plt
>>> plt.plot(H, fit)
>>> plt.xlabel('Field (T)')
>>> plt.ylabel('IRM fit')
>>> plt.show()
```

---

## `Langevin`

```python
Langevin(alpha)
```

Langevin function

**Parameters**

- **alpha** (`float`) — Langevin alpha value

**Returns**

- `float` — Langevin function value

---

## `Me_drift_correction`

```python
Me_drift_correction(H, M, descending_first=True)
```

Perform default IRM drift correction for a hysteresis loop based on the Me method.

This function applies a drift correction algorithm to magnetization data (M) measured as a function of applied field (H),
commonly used for IRM (Isothermal Remanent Magnetization) experiments. The correction is based on the Me signal,
which is the sum of the upper and reversed lower branches of the hysteresis loop.
The correction method adapts depending on whether significant drift is detected in the high-field region.

The drift estimate depends on measurement-time order, and the arrays are
expected in canonical order (descending upper branch first, as produced
by `grid_hyst_loop`). For a loop originally measured from negative
saturation, pass descending_first=False (detected from the raw field
values with `measured_descending_first`) so the correction is applied in
true time order rather than with the opposite time sense.

**Parameters**

- **H** (`numpy.ndarray`) — Array of magnetic field values, in canonical (descending-upper- branch-first) order.
- **M** (`numpy.ndarray`) — Array of measured magnetization values corresponding to `H`.
- **descending_first** (`bool`) — Whether the loop was originally measured with the descending branch first (default True). Use `measured_descending_first` on the raw field values to determine this for a gridded loop.

**Returns**

- `numpy.ndarray` — Corrected magnetization values after drift correction. 

**Examples**

```python
>>> H = np.linspace(-1, 1, 200)
>>> M = measure_hysteresis(H)
>>> M_cor = Me_drift_correction(H, M)
>>> plot(H, M, label='Original')
>>> plot(H, M_cor, label='Drift Corrected')
```

---

## `SD_MD_mixture`

```python
SD_MD_mixture(Mr_Ms_SD=0.5, Mr_Ms_MD=0.019, Bc_SD=400, Bc_MD=43, Bcr_SD=500, Bcr_MD=230, X_sd=0.6, X_MD=0.209, Xr_SD=0.48, Xr_MD=0.039)
```

function to calculate the SD/MD mixture curve according to Dunlop (2002)

**Parameters**

- **Mr_Ms_SD** (`float`) — remanent to saturation magnetization ratio for SD. The default is 0.5.
- **Mr_Ms_MD** (`float`) — remanent to saturation magnetization ratio for MD. The default is 0.019.
- **Bc_SD** (`float`) — coercivity for SD. The default is 400.
- **Bc_MD** (`float`) — coercivity for MD. The default is 43.
- **Bcr_SD** (`float`) — coercivity of remanence for SD. The default is 500.
- **Bcr_MD** (`float`) — coercivity of remanence for MD. The default is 230.
- **X_sd** (`float`) — approximate Mrs/Bc slope for SD. The default is 0.6.
- **X_MD** (`float`) — approximate Mrs/Bc slope for MD. The default is 0.209.
- **Xr_SD** (`float`) — approximate Mrs/Bcr slope for SD. The default is 0.48.
- **Xr_MD** (`float`) — approximate Mrs/Bcr slope for MD. The default is 0.039.

**Returns**

- `numpy.ndarray` — coercivity ratio array
- `numpy.ndarray` — saturation magnetization ratio array 
- `* the default values are fro the IRM database` — 

---

## `SP_SD_mixture`

```python
SP_SD_mixture(SP_size, SD_Mr_Ms=0.5, SD_Bcr_Bc=1.25, X_sd=3, T=300)
```

function to calculate the SP/SD mixture curve according to Dunlop (2002)

**Parameters**

- **SP_size** (`float`) — size of the superparamagnetic particle in nm
- **SD_Mr_Ms** (`float`) — remanent to saturation magnetization ratio. The default is 0.5.
- **SD_Bcr_Bc** (`float`) — remanent coercivity to coercivity ratio. The default is 1.25.
- **X_sd** (`float`) — approximate Mrs/Bc slope. The default is 3 for magnetite
- **T** (`float`) — temperature in Kelvin. The default is 300.

**Returns**

- `numpy.ndarray` — coercivity ratio array
- `numpy.ndarray` — saturation magnetization ratio array

---

## `SP_saturation_curve`

```python
SP_saturation_curve(SD_Mr_Ms=0.5, SD_Bcr_Bc=1.25)
```

function to calculate the SP saturation curve according to Dunlop (2002)

**Parameters**

- **SD_Mr_Ms** (`float`) — saturation magnetization ratio. The default is 0.5.
- **SD_Bcr_Bc** (`float`) — remanence coercivity to coercivity ratio. The default is 1.25.

**Returns**

- `numpy.ndarray` — coercivity ratio array
- `numpy.ndarray` — saturation magnetization ratio array

---

## `add_Bcr_to_specimens_table`

```python
add_Bcr_to_specimens_table(specimens_df, experiment_name, Bcr)
```

Add the Bcr value to the MagIC specimens table
the controled vocabulary for backfield derived Bcr is rem_bcr

**Parameters**

- **specimens_df** (`pandas.DataFrame`) — The specimens table from the MagIC database
- **experiment_name** (`str`) — The name of the experiment to which the Bcr value belongs
- **Bcr** (`float`) — The Bcr value to be added to the specimens table

---

## `add_curie_estimates_to_specimens_table`

```python
add_curie_estimates_to_specimens_table(specimens_df, experiment_name, estimates, method='inflection', branch='heating', critical_temp_type='Curie')
```

Write a Curie temperature estimate to a MagIC specimens table.

Sets ``critical_temp`` (in Kelvin, per the MagIC data model) and
``critical_temp_type`` (controlled vocabulary; 'Curie' by default) for
the rows whose ``experiments`` column matches ``experiment_name``, and
records the estimation method, branch, and uncertainty in the
``description`` column so the processing choice is archived with the
result. Updates ``specimens_df`` in place.

**Parameters**

- **specimens_df** (`pandas.DataFrame`) — MagIC specimens table (with an 'experiments' column).
- **experiment_name** (`str`) — Experiment the estimate belongs to.
- **estimates** (`pandas.DataFrame`) — Tidy table from ``curie_temperature_estimates``.
- **method** (`str`) — Which method's estimate to write (default 'inflection').
- **branch** (`str`) — Which branch's estimate to write (default 'heating').
- **critical_temp_type** (`str`) — MagIC controlled-vocabulary temperature type (default 'Curie').

---

## `add_hyst_stats_to_specimens_table`

```python
add_hyst_stats_to_specimens_table(specimens_df, hyst_results, overwrite=True)
```

Return a copy of the specimens table with hysteresis results added.

The input DataFrame is not modified. Assign the return value to
update your table, e.g.:
    specimens = add_hyst_stats_to_specimens_table(specimens, hyst_results)

**Parameters**

- **specimens_df** (`pandas.DataFrame`) — dataframe with the specimens data
- **hyst_results** (`pandas.DataFrame`) — DataFrame with hysteresis results including 'specimen' and 'experiment' columns, as output from rmag.process_hyst_loops. Has a numeric index (one row per experiment).
- **overwrite** (`bool`) — If True (default), existing MagIC column values and description stats are replaced with new values from hyst_results. If False, existing rows are preserved as-is and new rows are appended with the hyst results.

**Returns**

- `pandas.DataFrame` — A new DataFrame with hysteresis results added. If a specimen has multiple experiments, its row is duplicated so that each experiment gets its own row.

---

## `add_unmixing_to_specimens_table`

```python
add_unmixing_to_specimens_table(specimens_df, components_df, mode='rows')
```

Record coercivity unmixing results in a MagIC specimens table.

Two recording conventions are supported:

- mode='rows' (default, conformant with the MagIC data model): one new
  specimens row is appended per experiment and component, populating
  the controlled-vocabulary columns rem_cmf (component median field,
  in tesla), rem_cd (component dispersion, log10 units), and
  rem_n_comp (number of components in the model), along with specimen
  identity columns copied from the matching existing row. The full
  parameter set (proportion, skew, confidence intervals, ...) is
  stored as JSON in each new row's 'description' cell. Rows from a
  previous call for the same experiments are replaced.
- mode='description': the matching existing specimens rows are updated
  in place, storing the complete unmixing model as JSON under the
  'coercivity_unmixing' key in 'description' (any existing free text
  is preserved with a 'text | json' convention readable by
  parse_specimen_description).

**Parameters**

- **specimens_df** (`pandas.DataFrame`) — MagIC specimens table.
- **components_df** (`pandas.DataFrame`) — Tidy components table from unmix_backfield_experiments.
- **mode** (`str`) — 'rows' or 'description'.

**Returns**

- `pandas.DataFrame` — The updated specimens table. With mode='description' the input table is also modified in place; with mode='rows' a new table is returned (pandas cannot append rows in place).

---

## `aggregate_by_class`

```python
aggregate_by_class(components, boundaries_mT, class_names=None, coercivity_column='B_mean_mT', proportion_column='proportion', contribution_column='contribution', curve_factor=2.0, group_column='experiment', passthrough=('specimen', 'Bcr_mT', 'r_squared'))
```

Aggregate fitted unmixing components into coercivity classes.

Sums each component's remanence proportion (and, if available, its
contribution) into classes defined by one or more coercivity cut points,
per experiment. This is the robust way to quantify a mineral assemblage
from an unmixing fit: because it integrates the fitted distribution
within coercivity bands, the result is insensitive to how many
components the optimizer used or how it split a single mineral, so long
as the mineral populations are separated by the cut points. Typical use
is a magnetite/hematite split at a single boundary, but any number of
classes is supported.

**Parameters**

- **components** (`pandas.DataFrame`) — Tidy component table, e.g. from unmix_backfield_experiments (one row per experiment and component).
- **boundaries_mT** (`float or sequence of float`) — Coercivity cut point(s) in mT. A single value gives two classes; a sequence of k values gives k+1 classes. A component is placed in the class between the cut points that bracket its coercivity; a component exactly on a boundary goes to the lower class.
- **class_names** (`sequence of str`) — Names for the classes, in order of increasing coercivity (length = number of boundaries + 1). Defaults to 'class_1', 'class_2', ...; for a two-class split you would typically pass e.g. ['magnetite', 'hematite'].
- **coercivity_column** (`str`) — Column used to classify each component (default 'B_mean_mT').
- **proportion_column** (`str`) — Column summed to give each class's remanence fraction (default 'proportion').
- **contribution_column** (`str`) — Column summed (and divided by curve_factor) to give each class's absolute remanence; skipped if the column is absent.
- **curve_factor** (`float`) — Divisor applied to summed contributions. Use 2 for a shift-corrected backfield curve (which spans twice the saturation remanence) and 1 for an IRM acquisition curve.
- **group_column** (`str`) — Column identifying each specimen/experiment to aggregate within (default 'experiment').
- **passthrough** (`sequence of str`) — Columns copied through unchanged (first value per group), e.g. specimen name and fit statistics. Missing columns are ignored.

**Returns**

- `pandas.DataFrame` — One row per group with the passthrough columns and, for each class, a '{name}_fraction' column and (when contributions are present) a '{name}_remanence' column.

---

## `build_symmetric_hyst_grid`

```python
build_symmetric_hyst_grid(upper_branch, lower_branch)
```

Build a symmetric field grid over the overlap shared by the two loop branches.

---

## `calc_Bc`

```python
calc_Bc(H, M)
```

function for calculating the coercivity of the ferromagnetic component of a hysteresis loop
    the final Bc value is calculated as the average of the positive and negative Bc values

**Parameters**

- **H** (`numpy array`) — field values
- **M** (`numpy array`) — magnetization values

**Returns**

- `float` — coercivity of the ferromagnetic component of the hysteresis loop

---

## `calc_Mr_Mrh_Mih_Brh`

```python
calc_Mr_Mrh_Mih_Brh(grid_field, grid_magnetization)
```

function to calculate the Mrh and Mih values from a hysteresis loop

**Parameters**

- **grid_field** (`numpy array`) — gridded field values
- **grid_magnetization** (`numpy array`) — gridded magnetization values

**Returns**

- `numpy array` — field values of the upper branch (the two branches should have the same field values)
- `float` — remanent magnetization (Mrh interpolated at zero field)
- `numpy array` — remanent hysteretic magnetization, (upper - lower)/2
- `numpy array` — induced hysteretic magnetization, (upper + lower)/2
- `numpy array` — error curve err(H), the mismatch between the upper branch and the inverted lower branch
- `float` — median field of Mrh (field at which Mrh falls to half of Mr)

---

## `calc_Q`

```python
calc_Q(H, M, type='Q')
```

Calculate the quality factor (Q) for a magnetic hysteresis loop.

The Q factor is a logarithmic measure (base 10) of the signal-to-noise ratio for a
hysteresis loop, following Jackson and Solheid (2010): the upper and inverted lower
branches are treated as replicate measurements, so their mean squared moment relative
to the mean squared mismatch between them (the err(H) curve) quantifies signal/noise.
Q = log10(s/n); loops with Q >= 2 have small deviations from inversion symmetry while
loops with Q below ~0.3 (s/n ~ 2) are too noisy for meaningful parameter estimation.
The quality factor of the ferromagnetic component (Q_f of Jackson and Solheid, 2010)
is obtained by calling this function on the slope-corrected loop.

The calculation can be performed in two modes:
    - 'Q': Uses the mean squared magnetization of both the upper and lower branches.
    - 'Qf': Uses only the upper branch.

**Parameters**

- **H** (`array_like`) — Array of applied magnetic field values.
- **M** (`array_like`) — Array of measured magnetization (moment) values, corresponding to `H`.
- **type** (`(Q, Qf)`) — Type of Q calculation to perform:     - 'Q' (default): Uses both upper and lower branches of the loop.     - 'Qf': Uses only the upper branch.

**Returns**

- `float` — The calculated signal-to-noise ratio (before applying the logarithm).
- `float` — The quality factor, defined as log10(M_sn). 

**Notes**

- The function splits the hysteresis loop into upper and lower branches using `split_hyst_loop`.
- For type 'Q', the numerator is the average of the sum of squares of the upper and lower branches; for 'Qf', only the upper branch is used.
- The denominator is the sum of squares of err(H), the mismatch between the upper branch
  and the inverted lower branch, so M_sn is equivalent to the 1/(1 - R^2) signal/noise
  measure of Jackson and Solheid (2010, equation 3).
- Higher Q values indicate a higher signal-to-noise ratio in the hysteresis loop data.

**Examples**

```python
>>> H = np.linspace(-1, 1, 200)
>>> M = np.tanh(3 * H) + 0.05 * np.random.randn(200)
>>> M_sn, Q = calc_Q(H, M, type='Q')
>>> print(f"Signal-to-noise ratio: {M_sn:.3f}, Q: {Q:.2f}")
```

---

## `calc_verwey_estimate`

```python
calc_verwey_estimate(temps, mags, t_range_background_min=50, t_range_background_max=250, excluded_t_min=75, excluded_t_max=150, poly_deg=3)
```

Estimate the Verwey transition temperature and remanence loss of magnetite from MPMS data.
Plots the magnetization data, background fit, and resulting magnetite curve, and 
optionally the zero-crossing.

**Parameters**

- **temps** (`pd.Series`) — Series representing the temperatures at which magnetization measurements were taken.
- **mags** (`pd.Series`) — Series representing the magnetization measurements.
- **t_range_background_min** (`int or float`) — Minimum temperature for the background fitting range. Default is 50.
- **t_range_background_max** (`int or float`) — Maximum temperature for the background fitting range. Default is 250.
- **excluded_t_min** (`int or float`) — Minimum temperature to exclude from the background fitting range. Default is 75.
- **excluded_t_max** (`int or float`) — Maximum temperature to exclude from the background fitting range. Default is 150.
- **poly_deg** (`int`) — Degree of the polynomial for background fitting. Default is 3.

---

## `calc_zero_crossing`

```python
calc_zero_crossing(dM_dT_temps, dM_dT)
```

Calculate the temperature at which the second derivative of magnetization with respect to 
temperature crosses zero. This value provides an estimate of the peak of the derivative 
curve that is more precise than the maximum value.

The function computes the second derivative of magnetization (dM/dT) with respect to 
temperature, identifies the nearest points around the maximum value of the derivative, 
and then calculates the temperature at which this second derivative crosses zero using 
linear interpolation.

Parameters:
    dM_dT_temps (pd.Series): A pandas Series representing temperatures corresponding to
                             the first derivation of magnetization with respect to temperature.
    dM_dT (pd.Series): A pandas Series representing the first derivative of 
                       magnetization with respect to temperature.
Returns:
    float: The estimated temperature at which the second derivative of magnetization 
           with respect to temperature crosses zero.

Note:
    The function assumes that the input series `dM_dT_temps` and `dM_dT` are related to 
    each other and are of equal length.

---

## `chi_SP`

```python
chi_SP(SP_size, T)
```

SP size distribution function

**Parameters**

- **SP_size** (`float`) — size of the superparamagnetic particle in nm
- **T** (`float`) — temperature in Kelvin

**Returns**

- `float` — susceptibility value

---

## `clean_out_na`

```python
clean_out_na(dataframe)
```

Cleans a DataFrame by removing columns and rows that contain only NaN values.

Args:
    dataframe (pd.DataFrame): The DataFrame to be cleaned.
Returns:
    pd.DataFrame: A cleaned DataFrame with all-NaN columns and rows removed.

---

## `coercivity_curve_components`

```python
coercivity_curve_components(x, parameters, curve_type='backfield')
```

Evaluate each component in measurement space (cumulative curves).

For 'backfield' curves (processed so that magnetization decays from a
maximum toward zero with increasing field magnitude) each component is
contribution * (1 - CDF); for 'acquisition' curves each component is
contribution * CDF.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **parameters** (`pandas.DataFrame or array - like`) — Component parameters (see coercivity_spectrum_components).
- **curve_type** (`str`) — 'backfield' or 'acquisition'.

**Returns**

- `numpy.ndarray` — Array of shape (n_components, len(x)).

---

## `coercivity_curve_model`

```python
coercivity_curve_model(x, parameters, offset=0.0, curve_type='backfield')
```

Evaluate the summed unmixing model in measurement space.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **parameters** (`pandas.DataFrame or array - like`) — Component parameters (see coercivity_spectrum_components).
- **offset** (`float`) — Constant baseline added to the model (accounts for a small unsaturated or instrumental offset; default 0).
- **curve_type** (`str`) — 'backfield' or 'acquisition'.

**Returns**

- `numpy.ndarray` — Total model curve at x.

---

## `coercivity_prior_table`

```python
coercivity_prior_table(minerals=None)
```

Summarize the coercivity-component prior library as a table.

**Parameters**

- **minerals** (`list of str`) — Component names, each a key of COERCIVITY_COMPONENT_LIBRARY. Defaults to the whole library, ordered by central (geometric-mean) coercivity.

**Returns**

- `pandas.DataFrame` — One row per component with its family, mean-coercivity window (mT), dispersion window (dp, log10 units), skew (Azzalini alpha) window, the implied distribution asymmetry, and the leading literature source.

---

## `coercivity_spectrum_components`

```python
coercivity_spectrum_components(x, parameters)
```

Evaluate each unmixing component in spectrum space (dM/dlog10 B).

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **parameters** (`pandas.DataFrame or array - like`) — One row per component with columns 'contribution' (area under the component in magnetization units), 'location', 'dp', 'skew'.

**Returns**

- `numpy.ndarray` — Array of shape (n_components, len(x)).

---

## `coercivity_spectrum_from_curve`

```python
coercivity_spectrum_from_curve(x, magnetization, curve_type='backfield')
```

Compute a finite-difference coercivity spectrum from a remanence curve.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT), monotonically increasing.
- **magnetization** (`array - like`) — Magnetization values at x (shifted to positive for backfield data, e.g. the 'magn_mass_shift' column from process_backfield_data).
- **curve_type** (`str`) — 'backfield' (decaying curve, spectrum = -dM/dx) or 'acquisition' (growing curve, spectrum = dM/dx).

**Returns**

- `tuple` — (x_mid, spectrum) where x_mid are midpoints between successive x values and spectrum is the centered finite-difference derivative.

---

## `coercivity_spectrum_model`

```python
coercivity_spectrum_model(x, parameters)
```

Evaluate the summed unmixing model in spectrum space.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **parameters** (`pandas.DataFrame or array - like`) — Component parameters (see coercivity_spectrum_components).

**Returns**

- `numpy.ndarray` — Total model spectrum at x.

---

## `coercivity_unmixing_interactive`

```python
coercivity_unmixing_interactive(x, magnetization, n_components=2, method='spectrum', curve_type='backfield', vary_skew=True, figsize=(9, 5))
```

Interactive widget for choosing initial unmixing parameters visually.

Initial parameter choices strongly influence nonlinear unmixing fits.
This widget shows the coercivity spectrum with a live model built from
per-component sliders (peak field in mT on a log scale, proportion of
the total remanence, dispersion DP, and skew). Sliders are seeded from
automatic peak detection. Pressing "Fit" runs the chosen optimizer
(spectrum- or measurement-space) starting from the current slider
values and overlays the optimized model.

*Important*: run `%matplotlib widget` in the notebook first so the
figure updates live.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT), e.g. 'log_dc_field' from process_backfield_data.
- **magnetization** (`array - like`) — Remanence curve values at x (e.g. 'magn_mass_shift').
- **n_components** (`int`) — Number of components (default 2).
- **method** (`str`) — 'spectrum' or 'curve' -- the fitting approach used by the Fit button.
- **curve_type** (`str`) — 'backfield' or 'acquisition'.
- **vary_skew** (`bool`) — Include skew sliders and let the fit vary skew (default True).
- **figsize** (`tuple`) — Figure size.

**Returns**

- `dict` — A live handle with keys 'initial_parameters' (DataFrame updated as sliders move; pass to the unmixing functions) and 'result' (filled with the standardized result dictionary after Fit is pressed).

---

## `collapse_hyst_field_plateaus`

```python
collapse_hyst_field_plateaus(field, magnetization)
```

Average consecutive repeated field steps into a single point.

---

## `compare_unmixing_models`

```python
compare_unmixing_models(results)
```

Compare unmixing fits with different numbers of components.

Builds a comparison table with information criteria and sequential
F-tests. All results must be fits of the same data with the same
method ('spectrum' or 'curve'); AIC/BIC values are only meaningful
relative to one another under that condition. The F-test compares each
model to the previous (simpler) one; a small p-value indicates the
additional component produces a statistically significant improvement.
As emphasized by Maxbauer et al. (2016) and Egli (2003), statistical
significance alone should not decide the number of components --
independent knowledge of the likely magnetic mineralogy should inform
the choice.

**Parameters**

- **results** (`list of dict`) — Result dictionaries from unmix_coercivity_spectrum or unmix_backfield_curve, typically with increasing n_components.

**Returns**

- `pandas.DataFrame` — One row per model with n_components, n_params, rss, r_squared, aic, bic, delta_aic, delta_bic, F, and p_value columns.

---

## `convert_temperature`

```python
convert_temperature(temp_array, input_unit, output_unit)
```

Convert temperatures between Kelvin and Celsius.

**Parameters**

- **temp_array** (`array - like`) — Temperatures in `input_unit`.
- **input_unit** (`(K, C)`) — Unit of the input temperatures.
- **output_unit** (`(K, C)`) — Desired unit for output temperatures.

**Returns**

- `numpy.ndarray` — Temperatures converted to `output_unit`. 

**Raises**

- `ValueError` — If `input_unit` or `output_unit` is not 'K' or 'C'.

---

## `curie_Ms_squared_extrapolation`

```python
curie_Ms_squared_extrapolation(T, M, fit_range=None, exponent=2.0, min_points=5, min_m=None)
```

Mean-field (Moskowitz, 1981) extrapolation of Ms**exponent to zero.

Fits ``M**exponent`` linearly in temperature over ``fit_range`` and returns
the x-intercept (where ``M**exponent = 0``) as the Curie temperature. This
is the extrapolation method of Moskowitz (1981,
doi:10.1016/0012-821X(81)90028-5) for saturation-magnetization (Ms-T)
thermomagnetic curves whose signal terminates below Tc — e.g.
titanomaghemites that chemically invert on heating before Ms reaches zero,
so a direct Ms = 0 crossing is unavailable.

Near Tc a second-order (mean-field / Landau-Belov) treatment gives
Js proportional to (Tc - T)**(1/2), so Js**2 is linear in T and
extrapolates to zero at Tc (Moskowitz, 1981, eqs. 4-5). This is the
magnetization-space analog of the Curie-Weiss inverse-susceptibility
method: there 1/chi rises linearly to a zero at theta; here Ms**2 falls
linearly to a zero at Tc, so ``curie_temp = -intercept/slope`` with the
slope required to be negative. (Moskowitz normalizes by Js0 = Js(T0); the
normalization does not change the x-intercept and is omitted here.)

The mean-field form is valid for T/Tc > 0.8 (roughly the uppermost ~100 C
below Tc), so restrict ``fit_range`` to the steep descent below where the
curve flattens out or the sample alters; report the window. Following
Moskowitz, ``exponent=2`` (critical exponent beta = 1/2) is the
theoretically justified, best-conditioned choice and minimized his fit
error (+/-5 C) on standards (Fe3O4 572 C, Ni 360 C, CrO2 133 C). For a
mean-field curve exponent=2 is exact; any other exponent introduces
curvature into the fitted data and moves the extrapolated Tc off the true
value. Moskowitz reported that a smaller exponent shifted Tc lower for his
irreversible titanomaghemite samples; on ideal and mean-field curves a
smaller exponent instead biases Tc high, so the direction of the bias
depends on the true curve shape. For a multiphase sample only the Curie
temperature of the highest-Tc phase is recovered.

**Parameters**

- **T** (`array - like`) — Temperatures, ascending (Celsius or Kelvin; Tc is returned in the same unit as the input, since the x-intercept is unit-invariant).
- **M** (`array - like`) — Saturation-magnetization (Ms) values.
- **fit_range** (`tuple of (float, float)`) — Temperature interval for the linear fit of ``M**exponent``. If None, the upper 40 percent of the temperature span is used — a starting guess only; inspect the fit and set the window explicitly (below the flattening/alteration temperature) for reported values.
- **exponent** (`float`) — Power to which ``M`` is raised before the linear fit (default 2.0, i.e. beta = 1/2). Any other exponent introduces curvature into the fit and moves Tc off the mean-field value (Moskowitz, 1981).
- **min_points** (`int`) — Minimum number of usable points in the window (default 5; floored at 3, since the covariance fit needs more than two points).
- **min_m** (`float`) — Exclude points with magnetization below this value from the fit (same units as ``M``). Useful for screening out a flattened, weak, or altered tail.

**Returns**

- `dict` — ``curie_temp`` (Tc = -intercept/slope, in the input temperature unit), ``curie_temp_stderr`` (1-sigma from the fit covariance), ``params`` with ``slope``, ``intercept``, ``r_squared``, ``n_points``, ``fit_range``, ``exponent``, a ``note`` when the fit is not usable (too few points or a non-negative slope), and a ``warning`` when the fitted points are dominated by repeated (flattened/altered) values, and ``diagnostics`` with ``T``, ``m_pow`` (= ``M**exponent``), the fitted points, and the extrapolated line.

---

## `curie_derivative_estimates`

```python
curie_derivative_estimates(T, y, t_range=None, smooth_window=0)
```

Derivative-based Curie temperature estimates from one thermomagnetic branch.

Two estimates are returned:

* ``inflection_temp`` — the inflection point of the curve, located as the
  zero crossing of the second derivative between its extrema (refined by
  linear interpolation), with the minimum of the first derivative as
  fallback. For in-field magnetization curves M(T), Landau theory places
  the Curie temperature at this inflection point, independent of the
  applied field (Fabian et al., 2013, doi:10.1029/2012GC004440).
* ``max_curvature_temp`` — the maximum of the second derivative
  ("maximum curvature" of the concave part of the heating curve). This is
  the classical practical definition of Ade-Hall et al. (1965) that is
  implemented in many software packages (including the legacy
  ``ipmag.curie``). On M(T) curves it lies systematically above the
  inflection-point Tc (typically 10-15 degrees C) and shifts with the
  strength of the applied field (Fabian et al., 2013). It is searched on
  the concave shoulder above the steepest descent, consistent with this
  definition.

Non-finite temperature or magnetization values are dropped before
differentiation, and the steepest descent is located over the interior of
the branch (the one-sided derivatives at the first and last points are the
most noise-prone). Both estimates are anchored to that steepest-descent
point, so isolated noise or structure in the flat tails does not capture
them.

Differentiation amplifies noise, so even a smoothed signal can yield a
ragged first derivative whose global minimum is a noise spike within a
broad, flat-bottomed transition rather than the true steepest descent.
Set ``smooth_window`` to smooth the first and second derivatives on the
same temperature scale used to smooth the signal (as the legacy
``ipmag.curie`` does), which locates the estimate at the center of the
transition instead of an arbitrary spike; ``curie_temperature_estimates``
passes its smoothing window through automatically. For multi-phase curves,
additionally use ``t_range`` to isolate the transition of interest, and
check the estimate against ``first_derivative_min_temp`` and the
``diagnostics`` arrays.

**Parameters**

- **T** (`array - like`) — Temperatures, ascending (Celsius or Kelvin; the returned temperatures are in the same unit as the input).
- **y** (`array - like`) — Magnetization or susceptibility values.
- **t_range** (`tuple of (float, float)`) — Restrict the analysis to temperatures within ``(t_min, t_max)``. Useful for isolating one Curie transition in a multi-phase curve or excluding low-temperature structure (e.g., a Hopkinson peak).
- **smooth_window** (`float`) — Width, in the temperature units of ``T``, of a moving-average window applied to the first and second derivatives before their extrema are located (default 0, no derivative smoothing). Recommended for noisy data; a good choice is the window used to smooth the signal. The ``diagnostics`` arrays reflect the smoothed derivatives actually used.

**Returns**

- `dict` — ``inflection_temp``, ``max_curvature_temp``, ``first_derivative_min_temp`` (the discrete minimum of dy/dT), ``zero_crossing_temps`` (all interpolated zero crossings of the second derivative between its extrema), and ``diagnostics`` with the arrays ``T``, ``dy_dT``, and ``d2y_dT2``.

---

## `curie_inverse_susceptibility`

```python
curie_inverse_susceptibility(T, chi, fit_range=None, min_points=5, min_chi=None)
```

Curie-Weiss (inverse susceptibility) estimate of the ordering temperature.

Above the Curie temperature the susceptibility of the paramagnetic phase
follows the Curie-Weiss law chi = C / (T - theta), so 1/chi is linear in T
and extrapolates to zero at the paramagnetic Curie temperature theta. A
straight line is fit to 1/chi within ``fit_range`` and
``curie_temp = -intercept/slope`` is returned.

This is the recommended quantitative approach for low-field
susceptibility X(T) curves (Petrovsky & Kapicka, 2006,
doi:10.1029/2006JB004507). Two caveats apply: (1) theta is an upper bound
on the Curie temperature (theta >= Tc, with the difference depending on
the strength of magnetic interactions); (2) the estimate is sensitive to
the choice of fitting window — the fit must be restricted to temperatures
where the signal is fully paramagnetic, above the steep decrease, and
where the (holder-corrected) susceptibility is still resolved above the
instrument's measurement resolution. Well above the transition the
paramagnetic signal commonly becomes so weak that successive readings
repeat identical values and 1/chi shows flat plateaus; including
such quantized points flattens the fitted slope and biases theta low —
a theta below an independently estimated Tc (e.g., the inflection point)
is a red flag for this. The function detects repeated quantized values in
the fitting window and attaches a warning; use ``min_chi`` or tighten
``fit_range`` to exclude resolution-limited points. Report the fitting
window with the estimate.

**Parameters**

- **T** (`array - like`) — Temperatures, ascending (Celsius or Kelvin; theta is returned in the same unit as the input).
- **chi** (`array - like`) — Susceptibility values (holder-corrected).
- **fit_range** (`tuple of (float, float)`) — Temperature interval for the linear fit of 1/chi. If None, the upper 20 percent of the temperature span is used — a starting guess only; inspect the fit and set the window explicitly for reported values.
- **min_points** (`int`) — Minimum number of usable points in the window (default 5).
- **min_chi** (`float`) — Exclude points with susceptibility below this value from the fit (same units as ``chi``). Useful for screening out resolution-limited values in the high-temperature tail.

**Returns**

- `dict` — ``curie_temp`` (theta), ``curie_temp_stderr`` (1-sigma from the fit covariance), ``params`` with slope, intercept, ``r_squared``, ``curie_constant`` (1/slope), ``n_points``, ``fit_range``, and a ``warning`` when the fitted points are dominated by repeated (quantized) susceptibility values, and ``diagnostics`` with the 1/chi arrays and fitted line.

---

## `curie_inverse_susceptibility_interactive`

```python
curie_inverse_susceptibility_interactive(experiment, temperature_column='meas_temp', magnetic_column='susc_chi_mass', temp_unit='C', input_unit='K', smooth_window=0, remove_holder=True, branch='heating', initial_fit_range=None, figsize=(6, 6))
```

Interactive (Bokeh) Curie-Weiss fit to inverse susceptibility.

Displays 1/chi versus temperature with a two-point fit line whose
endpoints can be dragged (Bokeh ``PointDrawTool``); the extrapolated
temperature where 1/chi reaches zero — the paramagnetic Curie temperature
theta (Petrovsky & Kapicka, 2006, doi:10.1029/2006JB004507) — updates
live below the plot.

This tool is for exploration: it helps identify the temperature interval
over which 1/chi is linear (fully paramagnetic). For reported values, use
``curie_inverse_susceptibility`` with the ``fit_range`` identified here,
so that the fit is reproducible.

**Parameters**

- **experiment** (`pandas.DataFrame`) — MagIC-formatted experiment DataFrame.
- **temperature_column** — 
- **magnetic_column** — 
- **temp_unit** — 
- **input_unit** — 
- **smooth_window** — As in ``curie_temperature_estimates``.
- **remove_holder** — As in ``curie_temperature_estimates``.
- **branch** (`(heating, cooling)`) — Branch to display (default 'heating').
- **initial_fit_range** (`tuple of (float, float)`) — Initial temperatures (in ``temp_unit``) for the two fit-line endpoints. Defaults to 70 percent and 100 percent of the branch's temperature range — a starting position only, meant to be dragged.
- **figsize** (`tuple`) — (width, height) in inches; height sets the Bokeh plot height.

**Returns**

- `None` — The interactive Bokeh application is displayed as a side effect; no value is returned. For a reproducible, reportable estimate use ``curie_inverse_susceptibility`` with the ``fit_range`` identified here. 

**Raises**

- `ValueError` — If the requested ``branch`` is not present in the experiment, or if it has fewer than two usable (finite, positive-susceptibility) points.

---

## `curie_landau_fit`

```python
curie_landau_fit(T, M, fit_range=None, temp_unit='C', tc_bounds=None, n_grid=40)
```

Fit the in-field Landau equation of state to a magnetization curve M(T).

The model is ``M(T) = M0 * m(tau; h)`` where m is the physical root of the
Landau equation of state ``m**3 + tau*m = h`` with ``tau = (T - Tc)/Tc``
in absolute temperature (Fabian et al., 2013, doi:10.1029/2012GC004440,
eqs. 3-5). The reduced field h rounds the transition and produces the
field-induced tail above Tc, so the fit uses the full curve without an
ad-hoc baseline. In the ``h = 0`` limit the model is
``M = M0*sqrt(1 - T/Tc)`` below Tc — the mean-field form underlying the
extrapolation method of Moskowitz (1981,
doi:10.1016/0012-821X(81)90028-5).

The fit provides the physically grounded Curie temperature (at the
inflection point of the in-field curve), with a formal 1-sigma
uncertainty. When ``fit_range`` excludes the transition (all data below
Tc), the fit operates in Moskowitz-style extrapolation mode: this is the
only option for runs that end below Tc, but the reliability of the
extrapolated Tc decays rapidly with the distance between the highest
measured temperature and Tc, and the formal uncertainty then
underestimates the true (model-dependence dominated) uncertainty.

Caveats: the model describes a single ferromagnetic phase near its
transition; admixed paramagnetic signal, multiple phases, or alteration
during heating violate it — restrict ``fit_range`` to isolate one
transition. Applied to low-field susceptibility the model does not
describe the dominant susceptibility mechanisms (see Fabian et al., 2013)
and results should be treated as qualitative.

**Parameters**

- **T** (`array - like`) — Temperatures, ascending, in ``temp_unit``.
- **M** (`array - like`) — Magnetization values.
- **fit_range** (`tuple of (float, float)`) — Temperature interval (in ``temp_unit``) used in the fit. Default uses all points.
- **temp_unit** (`(C, K)`) — Unit of the input temperatures (default 'C'). The reduced temperature is always formed in Kelvin internally; results are returned in ``temp_unit``.
- **tc_bounds** (`tuple of (float, float)`) — Bounds for Tc in Kelvin (default: (min(T)+1 K, 2000 K)). Tighten for extrapolation fits when independent constraints exist.
- **n_grid** (`int`) — Number of Tc values in the coarse initialization grid (default 40).

**Returns**

- `dict` — ``curie_temp`` and ``curie_temp_stderr`` in ``temp_unit``, ``params`` with ``M0``, ``h``, ``rss``, ``n_points``, ``fit_range``, ``tc_bounds``, and ``extrapolation`` (True when the highest fitted temperature is below the fitted Tc), and ``diagnostics`` with the fitted data and a dense model curve for plotting.

---

## `curie_temperature_estimates`

```python
curie_temperature_estimates(experiment, methods=None, temperature_column='meas_temp', magnetic_column='susc_chi_mass', temp_unit='C', input_unit='K', smooth_window=0, remove_holder=True, branches=('heating', 'cooling'), method_kwargs=None, print_estimates=False, data_type=None, branch_data=None, return_method_results=False)
```

Estimate the Curie temperature of a thermomagnetic experiment with
multiple methods and return a tidy comparison table.

The experiment is preprocessed with ``prepare_thermomag_branches`` and
each requested method is applied to each requested branch. Systematic
differences between the estimates are expected and diagnostic: see the
module notes above and Lattard et al. (2006, doi:10.1029/2006JB004591)
for the magnitude of inter-method offsets on synthetic titanomagnetites.

The available methods are:

* ``'inflection'`` — inflection point (``curie_derivative_estimates``);
  the recommended estimator for in-field M(T).
* ``'max_curvature'`` — second-derivative maximum
  (``curie_derivative_estimates``); classical, biased high on M(T).
* ``'two_tangent'`` — intersecting tangents (``curie_two_tangent``);
  for M(T), discouraged on susceptibility.
* ``'inverse_susceptibility'`` — Curie-Weiss extrapolation
  (``curie_inverse_susceptibility``); for susceptibility, heating branch.
* ``'landau'`` — in-field Landau equation-of-state fit
  (``curie_landau_fit``); for M(T), supports extrapolation from runs
  that end below Tc.
* ``'ms_squared_extrapolation'`` — mean-field extrapolation of Ms^2 to
  zero (``curie_Ms_squared_extrapolation``; Moskowitz, 1981); for M(T)
  curves that end below Tc. Not selectable for susceptibility data.

**Parameters**

- **experiment** (`pandas.DataFrame`) — MagIC-formatted experiment DataFrame.
- **methods** (`sequence of str`) — Methods to apply. Default depends on the data type (see ``data_type``): susceptibility data use ``('inflection', 'max_curvature', 'inverse_susceptibility')``; magnetization data use ``('inflection', 'max_curvature', 'two_tangent', 'landau')``.
- **temperature_column** (`str`) — Name of the temperature column (default 'meas_temp').
- **magnetic_column** (`str`) — Name of the magnetization/susceptibility column (default 'susc_chi_mass').
- **temp_unit** (`(C, K)`) — Unit for reported temperatures (default 'C').
- **input_unit** (`(K, C)`) — Unit of the temperatures in ``experiment`` (default 'K', the MagIC convention).
- **smooth_window** (`float`) — Smoothing window width in ``temp_unit`` (default 0, no smoothing). Derivative-based methods generally require smoothing of noisy data; see ``optimize_moving_average_window``.
- **remove_holder** (`bool`) — Subtract the per-branch minimum (default True). Disable for runs that end below the Curie temperature.
- **branches** (`sequence of str`) — Branches to analyze, from ('heating', 'cooling').
- **method_kwargs** (`dict`) — Per-method keyword arguments, e.g. ``{'inverse_susceptibility': {'fit_range': (620, 700)}, 'landau': {'fit_range': (300, 650)}}``. The keys ``'inflection'`` and ``'max_curvature'`` (or the shared key ``'derivative'``) forward options such as ``t_range`` to ``curie_derivative_estimates``. Unknown keys raise a ValueError rather than being silently ignored.
- **print_estimates** (`bool`) — Print a one-line summary per estimate (default False).
- **data_type** (`(susceptibility, magnetization, None)`) — Explicit data type, controlling the default method set and the caveat notes. When None (default), inferred from whether ``magnetic_column`` contains 'susc' or 'chi'.
- **branch_data** (`dict`) — Precomputed output of ``prepare_thermomag_branches`` (with matching preprocessing arguments). When provided, preprocessing is skipped — used by ``plot_curie_estimates`` to avoid recomputation.
- **return_method_results** (`bool`) — If True, also return a dict keyed by ``(branch, method)`` holding each estimator's full result (including ``diagnostics``), for plotting or further analysis (default False).

**Returns**

- `pandas.DataFrame or (pandas.DataFrame, dict)` — One row per (branch, method) with columns ``specimen``, ``experiment``, ``branch``, ``method``, ``curie_temp``, ``curie_temp_stderr``, ``temp_unit``, ``params`` (dict of method-specific parameters), and ``notes``. With ``return_method_results=True``, additionally the per-(branch, method) result dicts.

---

## `curie_two_tangent`

```python
curie_two_tangent(T, y, lower_range=None, upper_range=None, min_points=3)
```

Two-tangent (intersecting tangents) Curie temperature estimate.

Straight lines are fit to a segment of the steeply descending limb below
the Curie temperature and to the near-flat baseline above it; the
temperature of their intersection is returned. The construction follows
Gromme et al. (1969, doi:10.1029/JB074i022p05277), who applied it to
strong-field J-T curves.

Applicability and bias: the method is intended for magnetization curves
M(T). Even there it coincides with the maximum-curvature estimate and
therefore lies systematically above the inflection-point Curie temperature
(Fabian et al., 2013, doi:10.1029/2012GC004440). Applied to low-field
susceptibility X(T) it lacks a rigorous physical basis and can overestimate
Tc (Petrovsky & Kapicka, 2006, doi:10.1029/2006JB004507); it remains a
robust, transition-based estimator useful for mineral identification and
for comparison with legacy results.

**Parameters**

- **T** (`array - like`) — Temperatures, ascending (Celsius or Kelvin; the returned intersection temperature is in the same unit as the input).
- **y** (`array - like`) — Magnetization (preferred) or susceptibility values.
- **lower_range** (`tuple of (float, float)`) — Temperature interval for the descending-limb tangent. If None, the contiguous region around the steepest descent where the slope is at least half the steepest slope is used.
- **upper_range** (`tuple of (float, float)`) — Temperature interval for the baseline tangent. If None, points above the steepest descent where the slope has decayed to within 5 percent of the steepest slope are used (falling back to the uppermost decile of points).
- **min_points** (`int`) — Minimum number of points required in each tangent segment (default 3).

**Returns**

- `dict` — ``curie_temp`` (intersection temperature, NaN if the tangents are parallel or a segment has too few points), ``params`` with the two (slope, intercept) pairs, the temperature ranges actually used, point counts, and ``reliable`` (False when no near-flat baseline was found above the transition and the upper tangent fell back to the uppermost points, i.e. the curve may end below Tc), and ``diagnostics`` with the segment masks for plotting.

---

## `estimate_coercivity_components`

```python
estimate_coercivity_components(x, spectrum, n_components, smooth_window=None)
```

Automatic initial-guess estimation for unmixing components.

The spectrum is interpolated onto a uniform grid, lightly smoothed
(Savitzky-Golay), and searched for peaks. The n_components most
prominent peaks seed the component locations; widths at half maximum
seed dp; peak heights seed the contributions. If fewer peaks than
components are found, the remaining components are placed at evenly
spaced quantiles of the cumulative spectrum.

Initial choices matter for nonlinear fitting: these automatic estimates
are a starting point that can (and often should) be refined by the user,
e.g. with coercivity_unmixing_interactive.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **spectrum** (`array - like`) — Coercivity spectrum values at x.
- **n_components** (`int`) — Number of components to estimate.
- **smooth_window** (`int`) — Savitzky-Golay window length (grid points). Defaults to ~1/10 of the grid (minimum 5).

**Returns**

- `pandas.DataFrame` — Initial parameters with columns 'contribution', 'location', 'dp', 'skew' (skew = 0), sorted by location.

---

## `estimate_measurement_noise`

```python
estimate_measurement_noise(x, magnetization, curve_type='backfield')
```

Robustly estimate the measurement noise of a remanence curve.

Suppresses the smooth signal and returns a robust standard deviation of
the residual scatter, without any smoothing or model fitting. This
provides the noise scale needed to judge when an unmixing model fits "to
within the measurement noise" (see select_n_components).

For each interior point the smooth signal is removed by comparing the
point to the value predicted for it by cubic interpolation through its
four nearest neighbours (two on each side) at their actual field
positions; the residual is scaled by the standard deviation that
combination would have under white noise, and a robust (median-absolute)
scale is taken. Because the interpolation is exact for any polynomial of
degree <= 3, this removes not just the local slope but the curvature of
the sigmoidal curve, so the estimator is little biased even on the coarse
and non-uniform (typically log-spaced) field grids of backfield/IRM
measurements -- unlike a plain second difference, which cancels only a
linear trend and inflates the noise where the curve bends.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT); the data are sorted by x internally.
- **magnetization** (`array - like`) — Remanence curve values at x.
- **curve_type** (`str`) — Unused; accepted for signature consistency with the unmixing functions.

**Returns**

- `float` — Estimated measurement noise standard deviation in magnetization units.

---

## `experiment_selection`

```python
experiment_selection(measurements, experiment_name)
```

This function filters a measurements DataFrame to return only the rows 
that correspond to the specified experiment name. 

**Parameters**

- **measurements** (`pd.DataFrame`) — The DataFrame containing measurement data with an 'experiment' column.
- **experiment_name** (`str`) — The name of the experiment to select from the DataFrame.

**Returns**

- `pd.DataFrame` — A DataFrame containing only the rows corresponding to the specified experiment.

---

## `extract_hyst_data`

```python
extract_hyst_data(df, specimen_name)
```

Extracts hysteresis loop data for a specific specimen from a dataframe.

This function filters measurements for a given specimen and returns rows
whose MagIC method codes contain 'LP-HYS' (hysteresis loop experiments).

Parameters:
    df (pandas.DataFrame): The dataframe containing MagIC measurement data.
    specimen_name (str): The name of the specimen to filter data for.

Returns:
    pandas.DataFrame: Measurements for the specimen with 'LP-HYS' in
        their method codes (empty if none are present).

Example:
    >>> hyst_data = extract_hyst_data(measurements_df, 'Specimen_1')

---

## `extract_mpms_data_dc`

```python
extract_mpms_data_dc(measurements, specimen_name)
```

Extracts and separates MPMS data for a specified specimen from a DataFrame.

This function filters data for the given specimen and separates it based on
MagIC measurement method codes. It specifically extracts data corresponding to
'LP-FC' (Field Cooled), 'LP-ZFC' (Zero Field Cooled), 'LP-CW-SIRM:LP-MC' (Room
Temperature SIRM measured upon cooling), and 'LP-CW-SIRM:LP-MW' (Room Temperature
SIRM measured upon Warming). For each method code, if the data is not available,
an empty DataFrame with the same columns as the specimen data is returned.

Parameters:
    measurements (pd.DataFrame): DataFrame containing MPMS measurement data.
    specimen_name (str): Name of the specimen to filter data for.

Returns:
    tuple: A tuple containing four DataFrames:
        - fc_data: Data filtered for 'LP-FC' method if available, otherwise an empty
          DataFrame.
        - zfc_data: Data filtered for 'LP-ZFC' method if available, otherwise an empty
          DataFrame.
        - rtsirm_cool_data: Data filtered for 'LP-CW-SIRM:LP-MC' method if available,
          otherwise an empty DataFrame.
        - rtsirm_warm_data: Data filtered for 'LP-CW-SIRM:LP-MW' method if available,
          otherwise an empty DataFrame.

---

## `find_hyst_turning_point`

```python
find_hyst_turning_point(field)
```

Find the single loop reversal, tolerating repeated plateaus and minor field glitches.

Returns the index of the last point of the first field sweep, such that
field[:index+1] is the first branch and field[index+1:] is the second
branch. Zero field steps (plateaus from repeated field values) are ignored
when detecting the reversal, and if field glitches produce multiple sign
changes, the reversal closest to the global field extremum opposite the
initial sweep direction is chosen.

---

## `goethite_removal`

```python
goethite_removal(rtsirm_warm_data, rtsirm_cool_data, t_min=150, t_max=290, poly_deg=2, rtsirm_cool_color='#17becf', rtsirm_warm_color='#d62728', symbol_size=4, return_data=False)
```

Analyzes and visualizes the removal of goethite signal from Room Temperature Saturation
Isothermal Remanent Magnetization (RTSIRM) warming and cooling data. The function fits
a polynomial to the RTSRIM warming curve between specified temperature bounds to model
the goethite contribution, then subtracts this fit from the original data. The corrected
and uncorrected magnetizations are plotted, along with their derivatives, to assess the
effect of goethite removal.

Parameters:
    rtsirm_warm_data (pd.DataFrame): DataFrame containing 'meas_temp' and 'magn_mass' columns
                                     for RTSIRM warming data.
    rtsirm_cool_data (pd.DataFrame): DataFrame containing 'meas_temp' and 'magn_mass' columns
                                     for RTSIRM cooling data.
    t_min (int, optional): Minimum temperature for polynomial fitting. Default is 150.
    t_max (int, optional): Maximum temperature for polynomial fitting. Default is 290.
    poly_deg (int, optional): Degree of the polynomial to fit. Default is 2.
    rtsirm_cool_color (str, optional): Color code for plotting cooling data. Default is '#17becf'.
    rtsirm_warm_color (str, optional): Color code for plotting warming data. Default is '#d62728'.
    symbol_size (int, optional): Size of the markers in the plots. Default is 4.
    return_data (bool, optional): If True, returns the corrected magnetization data for both
                                  warming and cooling. Default is False.

Returns:
    Tuple[pd.Series, pd.Series]: Only if return_data is True. Returns two pandas Series
                                 containing the corrected magnetization data for the warming
                                 and cooling sequences, respectively.

---

## `goethite_removal_interactive`

```python
goethite_removal_interactive(measurements, specimen)
```

Display an interactive widget for fitting and visualizing goethite removal from low temperature remanence data.

This function creates an interactive interface that allows the user to select a specimen and adjust parameters
(temperature range and polynomial degree) for fitting the goethite component in RTSIRM (Room Temperature Saturation Isothermal Remanent Magnetization) warming and cooling curves.
The user can visually explore the effect of these parameters on the fit and resulting goethite removal, with real-time updated plots.

**Parameters**

- **measurements** (`pandas.DataFrame`) — Low temperature remanence measurement data containing temperature and magnetization information for multiple specimens.
- **specimen** (`str or ipywidgets.Dropdown`) — Specimen to analyze, given either as a plain specimen name or as a selection widget (e.g. from ``specimen_selection_interactive``); for a widget the current ``.value`` is read when this function runs, so rerun the cell after changing the dropdown.

**Notes**

- Uses `ipywidgets` for interactive controls and `matplotlib` for plotting.
- The temperature range for the goethite fit and the polynomial degree of the fit can be adjusted via sliders.
- A reset button allows restoration of default parameter values.
- Supporting functions such as `extract_mpms_data_dc` and `goethite_removal` are required for this function to operate.
- This function is intended to be used in a Jupyter notebook or similar interactive environment.

**Returns**

- `None` — The function displays interactive widgets and plots for goethite removal but does not return a value. 

**Examples**

```python
>>> goethite_removal_interactive(measurements_df, specimen_dropdown)
>>> goethite_removal_interactive(measurements_df, 'A73-7-1350-4B-01a')
Displays interactive sliders and plots for fitting goethite removal to the selected specimen's data.
```

---

## `grid_hyst_loop`

```python
grid_hyst_loop(field, magnetization)
```

function to grid a hysteresis loop into a regular grid
    with grid intervals equal to the average field step size calculated from the data

**Parameters**

- **field** (`numpy array or list`) — hysteresis loop field values
- **magnetization** (`numpy array or list`) — hysteresis loop magnetization values

**Returns**

- `numpy array` — gridded field values
- `numpy array` — gridded magnetization values

---

## `hyst_HF_nonlinear_optimization`

```python
hyst_HF_nonlinear_optimization(H, M, HF_cutoff, fit_type, initial_guess=[1, 1, -0.1, -0.1], bounds=([0, 0, -np.inf, -np.inf], [np.inf, np.inf, 0, 0]))
```

Optimize a high-field nonlinear fit

**Parameters**

- **H** (`numpy.ndarray`) — Array of field values.
- **M** (`numpy.ndarray`) — Array of magnetization values.
- **HF_cutoff** (`float`) — Fraction of max(|H|) defining the lower bound of the high-field region.
- **fit_type** (`(IRM, Fabian, Fabian_fixed_beta)`) — Type of nonlinear model to fit.
- **initial_guess** (`list of float`) — Initial parameter guess for the optimizer. Defaults to [1, 1, -0.1, -0.1]: χ_HF = 1, Mₛ = 1, a₁ = –0.1, a₂ = –0.1 (or α, β for Fabian).
- **bounds** (`tuple of array-like`) — Lower and upper bounds for each parameter. Defaults to ([0, 0, -∞, -∞], [∞, ∞, 0, 0]): - Lower: χ_HF ≥ 0, Mₛ ≥ 0, a₁ ≥ –∞, a₂ ≥ –∞   - Upper: χ_HF ≤ ∞, Mₛ ≤ ∞, a₁ ≤ 0, a₂ ≤ 0   (for Fabian, α and β follow the same positions/limits).

**Returns**

- `dict` — Fit results with keys: - 'chi_HF', 'Ms', 'a_1', 'a_2' (for IRM) or   'chi_HF', 'Ms', 'alpha', 'beta' (for Fabian variants) - 'Fnl_lin': float, F statistic for the improvement of the nonlinear fit   over a linear fit (Jackson and Solheid, 2010, equation 21); values above   ~3-3.5 indicate a statistically significant improvement for the   4-parameter models (for 'Fabian_fixed_beta' the degrees of freedom are   (1, N-3) and the 5% critical value is ~3.9-4.0). Because the nonlinear   coefficients are constrained non-positive, the statistic is conservative   under the null (saturated loops give values well below the critical   value, occasionally marginally negative when the bounded fit is a hair   worse than unconstrained least squares).

---

## `hyst_linearity_test`

```python
hyst_linearity_test(grid_field, grid_magnetization)
```

function for testing the linearity of a hysteresis loop

**Parameters**

- **grid_field** (`numpy array`) — gridded field values
- **grid_magnetization** (`numpy array`) — gridded magnetization values

**Returns**

- `dict` — dictionary of the results of the linearity test and intermediate statistics for the ANOVA calculation

---

## `hyst_loop_centering`

```python
hyst_loop_centering(grid_field, grid_magnetization)
```

function for finding the optimum applied field offset value for minimizing a linear fit through 
    the Me based on the R2 value. The idea is maximizing the residual noise in the Me gives the best centered loop. 

**Parameters**

- **grid_field** (`numpy array`) — gridded field values
- **grid_magnetization** (`numpy array`) — gridded magnetization values

**Returns**

- `float` — optimized applied field offset value for the loop
- `float` — calculated magnetization offset value for the loop based on the optimized applied field offset (intercept of the fitted line using the upper branch and the inverted and optimally offsetted lower branch)
- `float` — R-squared value of the linear fit between the upper branch and the inverted and offsetted lower branch

---

## `hyst_loop_centering_iterative`

```python
hyst_loop_centering_iterative(grid_field, grid_magnetization, hf_cutoff=0.8, low_field_fraction=0.35, shift_bound_fraction=0.1, weight_power=4, max_iterations=5, field_tolerance=1e-05, moment_tolerance=1e-08)
```

Center a hysteresis loop by iterating between provisional slope removal and offset fitting.

The routine alternates between fitting a provisional high-field slope and optimizing horizontal
and vertical offsets on the residual loop using a low-field-weighted inversion-symmetry metric.
This is designed for weak ferromagnetic loops superimposed on a strong linear background.

---

## `hyst_loop_saturation_test`

```python
hyst_loop_saturation_test(grid_field, grid_magnetization, max_field_cutoff=0.97)
```

Assess the saturation state of a magnetic hysteresis loop based on linearity at high-field segments.

This function evaluates the degree of saturation in a hysteresis loop by calculating the F statistic
for nonlinearity (FNL, the lack-of-fit F ratio of Jackson and Solheid, 2010) over high-field windows
starting at 60%, 70%, and 80% of the maximum field (up to a specified cutoff). A significant FNL
(above the 2.5 threshold) indicates reproducible curvature in that window, i.e. the ferromagnetic
moment has not saturated and a linear high-field fit is inappropriate there.

**Parameters**

- **grid_field** (`array_like`) — Array of applied magnetic field values for the hysteresis loop.
- **grid_magnetization** (`array_like`) — Array of magnetization (moment) values corresponding to `grid_field`.
- **max_field_cutoff** (`float`) — Fraction of the maximum field to use as an upper cutoff for the analysis (default is 0.97).

**Returns**

- `dict` — Dictionary containing:     - 'FNL60': float, FNL for the window from 60% of the maximum field.     - 'FNL70': float, FNL for the window from 70% of the maximum field.     - 'FNL80': float, FNL for the window from 80% of the maximum field.     - 'saturation_cutoff': float, lowest field fraction (0.6, 0.7, or 0.8) at which the       high-field segment is statistically linear (saturated); 0.92 (the IRM default for       a nonlinear fit window) if no tested window is linear.     - 'loop_is_saturated': bool, True if the loop is saturated (linear) in at least one       tested high-field window; False if all windows show significant nonlinearity,       in which case an approach-to-saturation fit should be used. 

**Notes**

- The function uses `loop_saturation_stats` to compute FNL values for each field fraction.
- FNL values below 2.5 indicate statistically linear (saturated) high-field behavior;
  values above 2.5 indicate significant nonlinearity (nonsaturation).
- The result is converted to standard Python types using `_to_native_python`.

**Examples**

```python
>>> results = hyst_loop_saturation_test(fields, magnetizations)
>>> print(results['saturation_cutoff'], results['loop_is_saturated'])
0.8 False
```

---

## `hyst_slope_correction`

```python
hyst_slope_correction(grid_field, grid_magnetization, chi_HF)
```

function for subtracting the paramagnetic/diamagnetic slope from a hysteresis loop
     the input should be gridded field and magnetization values

**Parameters**

- **grid_field** (`numpy array`) — gridded field values
- **grid_magnetization** (`numpy array`) — gridded magnetization values
- **chi_HF** (`float`) — X_HF

**Returns**

- `numpy array` — corrected ferromagnetic component of the magnetization

---

## `landau_magnetization`

```python
landau_magnetization(tau, h)
```

Reduced magnetization from the Landau equation of state with field term.

Solves ``m**3 + tau*m = h`` for the physical (largest real) root, where
``tau = (T - Tc)/Tc`` is the reduced temperature and ``h`` is the reduced
field (Fabian et al., 2013, doi:10.1029/2012GC004440, eqs. 3-5). For
``h = 0`` this reduces to ``m = sqrt(-tau)`` below the Curie temperature
and ``m = 0`` above it; for ``h > 0`` the transition is rounded and a
field-induced tail persists above Tc.

**Parameters**

- **tau** (`array - like`) — Reduced temperatures (T - Tc)/Tc; requires absolute temperatures (Kelvin) in the ratio.
- **h** (`float`) — Reduced field, >= 0.

**Returns**

- `numpy.ndarray` — Reduced magnetization m at each tau.

---

## `linear_HF_fit`

```python
linear_HF_fit(field, magnetization, HF_cutoff=0.8)
```

function to fit a linear function to the high field portion of a hysteresis loop

**Parameters**

- **field** (`numpy array or list`) — raw hysteresis loop field values
- **magnetization** (`numpy array or list`) — raw hysteresis loop magnetization values

**Returns**

- `float` — high-field susceptibility of the paramagnetic/diamagnetic contribution in SI units (the raw fitted slope in field units of Tesla multiplied by mu_0 = 4*pi*1e-7); `hyst_slope_correction` performs the inverse conversion when removing this contribution from a loop
- `float` — y-intercept of the linear fit can be interpreted to be the saturation magnetization of the ferromagnetic component

---

## `loop_closure_test`

```python
loop_closure_test(H, Mrh, HF_cutoff=0.8, Me=None, max_field_cutoff=0.99)
```

function for testing whether a hysteresis loop is closed at high fields

Mrh should be an even function of field for a well-behaved loop
(Mrh(-H) = Mrh(H)), so its field-reflection average (even part, with
unphysical negative values set to zero) is taken as the signal. A loop
that remains open at high fields (e.g. due to unsaturated high-coercivity
phases such as hematite or goethite) retains a significant Mrh signal in
the high-field window, giving a high signal-to-noise ratio (SNR) and a
high ratio of high-field Mrh area to total Mrh area (HAR). The noise is
estimated from the high-field portion of the err(H) curve when `Me` is
provided (matching the HystLab implementation of this test; Paterson et
al., 2018, section 4.5), or from the odd part of Mrh otherwise. Fields
above max_field_cutoff (default 99%) of the maximum field are excluded
from the high-field windows to avoid extreme-tip artifacts.

**Parameters**

- **H** — field values of the upper branch (ascending)
- **Mrh** — remanent hysteretic magnetization Mrh(H)
- **HF_cutoff** — high field cutoff value taken as fraction of the max field value
- **Me** — error curve err(H) on the same field axis (as returned by calc_Mr_Mrh_Mih_Brh); used as the noise estimate when provided
- **max_field_cutoff** — upper trim of the high-field windows as fraction of the max field

**Returns**

- `dict` — Dictionary containing:     - 'SNR': float, high-field signal-to-noise ratio in dB     - 'HAR': float, high-field to total Mrh area ratio in dB     - 'loop_is_closed': bool, True if SNR < 8 dB or HAR < -48 dB

---

## `loop_saturation_stats`

```python
loop_saturation_stats(field, magnetization, HF_cutoff=0.8, max_field_cutoff=0.97)
```

ANOVA statistics for the high field portion of a hysteresis loop

**Parameters**

- **field** (`numpy array`) — field values
- **magnetization** (`numpy array`) — magnetization values
- **HF_cutoff** (`float`) — high field cutoff value default is 0.8

**Returns**

- `dict` — dictionary of the results of the ANOVA calculation and intermediate statistics for the ANOVA calculation

---

## `magnetite_Ms`

```python
magnetite_Ms(T)
```

Magnetite saturation magnetization calculation

**Parameters**

- **T** (`float`) — temperature in Celsius

**Returns**

- `float` — saturation magnetization value

---

## `make_experiment_df`

```python
make_experiment_df(measurements, exclude_method_codes=None)
```

Creates a DataFrame of unique experiments from the measurements DataFrame.

**Parameters**

- **measurements** (`pd.DataFrame`) — The DataFrame containing measurement data with columns 'specimen',  'method_codes', and 'experiment'.
- **exclude_method_codes** (`list of str`) — List of method codes to exclude from the output DataFrame. Rows with  'method_codes' containing any of these substrings will be removed.

**Returns**

- `pd.DataFrame` — A DataFrame containing unique combinations of 'specimen', 'method_codes',  and 'experiment'.

---

## `measured_descending_first`

```python
measured_descending_first(field)
```

Whether the first-measured branch sweeps downward from positive field.

Gridding canonicalizes a loop to descending-upper-branch-first array
order regardless of how it was measured, so the original sweep order must
be detected from the raw field values before gridding and passed to the
time-order-sensitive drift corrections.

**Parameters**

- **field** (`array_like`) — Raw (ungridded) applied field values in measurement order.

**Returns**

- `bool` — True if the sweep starts at the positive field extreme (descending branch measured first), False if it starts at the negative extreme.

---

## `mineral_priors`

```python
mineral_priors(mineral_names, tighten=1.0, widen=1.0, field_max_mT=None, overrides=None)
```

Build a Bayesian-unmixing prior dictionary from named mineral components.

Converts a list of component names from COERCIVITY_COMPONENT_LIBRARY
into the `priors` dictionary accepted by unmix_coercivity_bayes, with a
mean-coercivity, dispersion, and skew window per component. The mean
window constrains the component's mean coercivity directly (the
skew-normal location is derived from the sampled dp and skew), so the
window is meaningful whatever the fitted skew. Components are sorted by
their central coercivity so the returned order matches the (low-to-high
coercivity) component ordering of the fit.

The library windows are broad starting points; real samples often need
them adapted. Use `tighten` to narrow every window, `widen` to broaden
every window, and `overrides` to replace specific windows outright when
a mineral in your samples sits outside its library range (for example a
harder, finer-grained magnetite, or a hematite whose low-coercivity
shoulder starts below the library's detrital-hematite window). A typical
workflow is to fit unconstrained first, see where the components land,
and then set the windows accordingly.

Because these are informative priors on an ill-posed decomposition, they
should be treated as soft constraints; see the caveats in the
COERCIVITY_COMPONENT_LIBRARY documentation.

**Parameters**

- **mineral_names** (`list of str`) — Component names, each a key of COERCIVITY_COMPONENT_LIBRARY (e.g. ['magnetite_detrital', 'hematite_pigmentary', 'hematite_detrital']).
- **tighten** (`float`) — Factor (>= 1) by which to narrow every window about its center; 1 keeps the library width.
- **widen** (`float`) — Factor (>= 1) by which to broaden every window about its center; 1 keeps the library width. tighten and widen compose (net width factor = widen / tighten).
- **field_max_mT** (`float`) — Maximum applied field of the measurement (mT). Coercivity upper bounds are clipped to this value, so that components whose window extends past the measured field (e.g. hard hematite or pyrrhotite in a 1-2 T experiment) are not given prior support the data cannot constrain.
- **overrides** (`dict`) — Per-mineral window replacements, mapping a mineral name to a dict with any of 'B_median_mT', 'dp', 'skew' as (low, high) tuples. The replacement window is used in place of the library value (before tighten/widen and field clipping are applied), letting you anchor to the library for most minerals while tuning the ones your samples require.

**Returns**

- `dict` — A priors dictionary with 'mean', 'dp', and 'skew' keys, each a list of (low, high) tuples in the order of increasing coercivity, suitable for `unmix_coercivity_bayes(..., priors=...)`. The 'mean' window constrains each component's MEAN coercivity (log10 mT) -- so a library window means what it says regardless of the component's skew. The chosen component names, in the same order, are returned under the 'components' key.

---

## `mpms_signal_blender`

```python
mpms_signal_blender(measurement_1, measurement_2, spec_1, spec_2, experiments=['LP-ZFC', 'LP-FC', 'LP-CW-SIRM:LP-MC', 'LP-CW-SIRM:LP-MW'], temp_col='meas_temp', moment_col='magn_mass', fraction=0.5)
```

function for simulating simple mixtures of MPMS dc remanence measurements using the Insitute for Rock Magnetism's
 rock magnetism bestiary data

**Parameters**

- **measurement_1** (`pandas.DataFrame`) — MagIC formatted dataframe containing the first set of measurements.
- **measurement_2** (`pandas.DataFrame`) — MagIC formatted dataframe containing the second set of measurements.
- **spec_1** (`str`) — Specimen name for the first set of measurements.
- **spec_2** (`str`) — Specimen name for the second set of measurements.
- **experiments** (`list of str`) — List of experiment method codes to consider for blending. Default is ['LP-ZFC', 'LP-FC', 'LP-CW-SIRM:LP-MC', 'LP-CW-SIRM:LP-MW'].
- **temp_col** (`str`) — Column name for temperature in the measurement dataframes. Default is 'meas_temp'.
- **moment_col** (`str`) — Column name for magnetization in the measurement dataframes. Default is 'magn_mass'.
- **fraction** (`float`) — Fraction of the first specimen's magnetization to blend with the second specimen's magnetization. Default is 0.5.

**Returns**

- `dict` — A dictionary where keys are experiment method codes and values are dictionaries containing:

---

## `mpms_signal_blender_interactive`

```python
mpms_signal_blender_interactive(measurement_1, measurement_2, experiments=['LP-ZFC', 'LP-FC', 'LP-CW-SIRM:LP-MC', 'LP-CW-SIRM:LP-MW'], temp_col='meas_temp', moment_col='magn_mass', figsize=(12, 6))
```

function for making interactive blender of MPMS dc remanence measurements using the Institute for Rock Magnetism's
 rock magnetism bestiary data

**Parameters**

- **measurement_1** (`pandas.DataFrame`) — MagIC formatted dataframe containing the first set of measurements.
- **measurement_2** (`pandas.DataFrame`) — MagIC formatted dataframe containing the second set of measurements.
- **experiments** (`list of str`) — List of experiment method codes to consider for blending. Default is ['LP-ZFC', 'LP-FC', 'LP-CW-SIRM:LP-MC', 'LP-CW-SIRM:LP-MW'].
- **temp_col** (`str`) — Column name for temperature in the measurement dataframes. Default is 'meas_temp'.
- **moment_col** (`str`) — Column name for magnetization in the measurement dataframes. Default is 'magn_mass'.
- **figsize** (`tuple of float`) — Size of the figure for plotting. Default is (12, 6).

---

## `optimize_moving_average_window`

```python
optimize_moving_average_window(experiment, min_temp_window=0, max_temp_window=50, steps=50, colormapwarm='tab20b', colormapcool='tab20c')
```

Visualize and optimize the moving average window size for smoothing experimental temperature-dependent data.

This function evaluates the effect of different moving average window sizes on the smoothing of both the warm and cool cycles
of an experiment (such as low temperature remanence or thermal demagnetization data). It iterates over a range of window sizes,
applies smoothing, and computes the average variance and root mean square (RMS) for each window. These metrics are plotted
to help the user visually identify the optimal window size for minimizing variance and RMS, balancing noise reduction and signal fidelity.

**Parameters**

- **experiment** (`object or structured array`) — Experimental data containing temperature and measurement values. It must be compatible with the `split_heating_cooling` function.
- **min_temp_window** (`float`) — Minimum window size (in degrees Celsius) for the moving average. Default is 0.
- **max_temp_window** (`float`) — Maximum window size (in degrees Celsius) for the moving average. Default is 50.
- **steps** (`int`) — Number of window size steps to evaluate between the minimum and maximum. Default is 50.
- **colormapwarm** (`str`) — Matplotlib colormap name for the warm cycle plot. Default is 'tab20b'.
- **colormapcool** (`str`) — Matplotlib colormap name for the cool cycle plot. Default is 'tab20c'.

**Returns**

- `matplotlib.figure.Figure` — The matplotlib Figure object containing the optimization plots.
- `numpy.ndarray of matplotlib.axes.Axes` — Array of Axes objects (one for the warm cycle, one for the cool cycle). 

**Examples**

```python
>>> fig, axs = optimize_moving_average_window(my_experiment, min_temp_window=5, max_temp_window=30, steps=20)
>>> fig.show()
```

---

## `parse_specimen_description`

```python
parse_specimen_description(description)
```

Parse a MagIC specimens 'description' cell into (text, dict).

The description convention used here stores free text and a
machine-readable JSON dictionary separated by ' | '. Legacy cells that
contain a Python dict repr (from older rockmagpy versions) are parsed
with ast.literal_eval.

**Parameters**

- **description** (`str or NaN`) — Contents of the description cell.

**Returns**

- `tuple` — (text, data) where text is the free-text portion (str, possibly empty) and data is the parsed dictionary (possibly empty).

---

## `plot_M_T`

```python
plot_M_T(data, temperature_column='meas_temp', magnetization_column='magn_mass', input_unit='K', plot_unit='K', interactive=False, return_figure=False, show_plot=True, size=(6, 3), legend_location='upper left')
```

Plot magnetization versus temperature in static or interactive mode.

**Parameters**

- **data** (`pandas.DataFrame or array - like`) — Table or array containing temperature and magnetization data.
- **temperature_column** (`str`) — Name of the temperature column in `data`.
- **magnetization_column** (`str`) — Name of the magnetization column in `data`.
- **input_unit** (`(K, C)`) — Unit of the input temperature data.
- **plot_unit** (`(K, C)`) — Unit for the x-axis display.
- **interactive** (`bool`) — If True, use Bokeh for an interactive plot.
- **return_figure** (`bool`) — If True, return the figure object(s). Assign to capture, e.g.:     fig, ax = plot_M_T(..., return_figure=True)
- **show_plot** (`bool`) — If True, display the plot immediately.
- **size** (`tuple(float, float)`) — Figure size in inches (Matplotlib) or height for Bokeh.
- **legend_location** (`str`) — Legend location in Matplotlib terms.

**Returns**

- `tuple or layout or None` — - If `return_figure=True` and `interactive=False`, returns (fig, ax). - If `return_figure=True` and `interactive=True`, returns the Bokeh   layout object. - Otherwise, returns None.

---

## `plot_backfield_data`

```python
plot_backfield_data(experiment, field='treat_dc_field', magnetization='magn_mass', Bcr=None, figsize=(5, 10), plot_raw=True, plot_processed=True, plot_spectrum=True, interactive=False, return_figure=False, show_plot=True, y_axis_units='Am²/kg', legend_location='upper left')
```

Plot backfield data: raw, processed, and coercivity spectrum.

**Parameters**

- **experiment** (`DataFrame`) — Must contain raw and, if requested, processed columns.
- **field** (`str`) — Name of the magnetic field column.
- **magnetization** (`str`) — Name of the magnetization column.
- **Bcr** (`float`) — Calculated Bcr (T). If provided, will be plotted as a pink star.
- **figsize** (`tuple(float, float)`) — Figure size (in inches).
- **plot_raw** (`bool`) — 
- **plot_processed** (`bool`) — 
- **plot_spectrum** (`bool`) — 
- **interactive** (`bool`) — 
- **return_figure** (`bool`) — 
- **show_plot** (`bool`) — 
- **y_axis_units** (`str`) — Units to display on the y-axis labels of raw and processed panels.
- **legend_location** (`str`) — Location of the legend in Matplotlib terms.

**Returns**

- `Matplotlib (fig, axes) or Bokeh grid or None` — 

---

## `plot_chi_T`

```python
plot_chi_T(experiment, temperature_column='meas_temp', magnetic_column='susc_chi_mass', temp_unit='C', smooth_window=0, remove_holder=True, plot_derivative=True, plot_inverse=False, interactive=True, return_figure=False, figsize=(6, 6), window_type='hanning')
```

Plot the high-temperature susceptibility curve, and optionally its derivative
and reciprocal using Bokeh or Matplotlib.

Parameters:
    experiment (pandas.DataFrame): MagIC-formatted experiment DataFrame.
    temperature_column (str): Name of temperature column.
    magnetic_column (str): Name of susceptibility column.
    temp_unit (str): "C" or "K" for the plotted temperatures (input
        temperatures are assumed to be in Kelvin, the MagIC convention).
    smooth_window (int): Window for smoothing, if 0, no smoothing is applied.
    remove_holder (bool): Subtract holder signal.
    plot_derivative (bool): Plot derivative.
    plot_inverse (bool): Plot inverse.
    interactive (bool): True for Bokeh, False for Matplotlib.
    return_figure (bool): Return figure objects if True.
    figsize (tuple): (width, height) in inches.
    window_type (str): Weighting function applied within each smoothing
        window, one of 'flat', 'hanning', 'hamming', 'bartlett', or
        'blackman' (default 'hanning'). Only used when smooth_window > 0.
        See prepare_thermomag_branches for a description of each option.

Returns:
    tuple or None: If return_figure is True, a tuple of the figure
    objects created (Matplotlib Figures when interactive is False, Bokeh
    figures when interactive is True). Otherwise the figures are displayed
    and None is returned.

---

## `plot_coercivity_prior_library`

```python
plot_coercivity_prior_library(minerals=None, field_range=(1.0, 5000.0), n_grid=400, figsize=None, ax=None)
```

Visualize the coercivity-component prior library.

Draws, for each named mineral component, the skew-normal coercivity
distribution implied by the centre of its library windows (mean
coercivity, dispersion, and skew) as a ridgeline over a shared log field
axis, with the mean-coercivity window drawn as a bar at the baseline and
a dot at its centre. Components are coloured by mineral family and ordered
by central coercivity, so the coercivity ranges of the different minerals,
their characteristic widths and skews, and -- importantly -- their
overlaps (the reason a coercivity prior constrains a window rather than
identifying a mineral) are all legible at a glance.

**Parameters**

- **minerals** (`list of str`) — Component names to show (default: the whole library).
- **field_range** (`tuple`) — (min, max) field in mT for the coercivity axis (default 1-5000).
- **n_grid** (`int`) — Number of points for the density curves.
- **figsize** (`tuple`) — Figure size; a default is chosen from the number of components.
- **ax** (`matplotlib.axes.Axes`) — Axis to draw on; a new figure is created if omitted.

**Returns**

- `tuple` — (fig, ax).

---

## `plot_coercivity_unmixing`

```python
plot_coercivity_unmixing(result, show_components=True, show_bootstrap=True, show_initial=False, n_grid=300, figsize=None, title=None, color_by='component', class_boundaries=None, class_colors=None)
```

Plot an unmixing result in its fitted data space.

For spectrum-space fits a single panel shows the data, total model, and
components. For measurement-space (curve) fits an upper panel shows the
measured curve with the cumulative model and a lower panel shows the
finite-difference spectrum of the data with the implied component
density curves. Bootstrap 95% bands are drawn when present.

**Parameters**

- **result** (`dict`) — Result from unmix_coercivity_spectrum, unmix_backfield_curve, or unmixing_bootstrap.
- **show_components** (`bool`) — Draw the individual components (default True).
- **show_bootstrap** (`bool`) — Draw bootstrap confidence bands if available (default True).
- **show_initial** (`bool`) — Also draw the model implied by the initial parameters (dotted), useful for judging how far the optimizer moved (default False).
- **n_grid** (`int`) — Number of points for smooth model curves.
- **figsize** (`tuple`) — Figure size; defaults depend on the number of panels.
- **title** (`str`) — Figure title.
- **color_by** (`str`) — How to color the components. 'component' (default) colors by fit order (C0, C1, ...). 'class' (or the alias 'coercivity') colors each component by the coercivity class its mean field falls into, so the mineralogy is read directly and a second component of the same mineral does not take a different color; requires class_boundaries.
- **class_boundaries** (`float or sequence of float`) — Coercivity cut points in mT that partition the components into classes when color_by='class' (e.g. 200 for a magnetite/hematite split, or [30, 300] for three classes).
- **class_colors** (`sequence`) — One color per class (length = number of boundaries + 1). Defaults to blue/red for two classes, or a diverging colormap otherwise.

**Returns**

- `tuple` — (fig, axes) with axes a list of the panel axes.

---

## `plot_curie_estimates`

```python
plot_curie_estimates(experiment, methods=None, temperature_column='meas_temp', magnetic_column='susc_chi_mass', temp_unit='C', input_unit='K', smooth_window=0, remove_holder=True, branches=('heating', 'cooling'), method_kwargs=None, figsize=(10, 10), return_figure=False, save_path=None, data_type=None)
```

Plot a thermomagnetic curve with Curie temperature estimates from
multiple methods overlain.

Produces a static matplotlib figure with up to three stacked panels:
(a) the (holder-corrected) curve per branch with a vertical line at each
method's estimate, the two-tangent construction, and the Landau model
curve where those methods are requested; (b) the first and second
derivatives with the inflection point and curvature maximum marked; and
(c) 1/chi with the Curie-Weiss fit when ``'inverse_susceptibility'`` is
requested. Method colors follow the colorblind-safe palette of Okabe &
Ito (2008); heating estimates are drawn with solid lines and cooling
estimates with dashed lines.

Parameters mirror ``curie_temperature_estimates``; see that function for
the estimation details and method caveats.

**Parameters**

- **experiment** (`pandas.DataFrame`) — MagIC-formatted experiment DataFrame.
- **methods** (`sequence of str`) — Methods to display (default as in ``curie_temperature_estimates``).
- **temperature_column** — 
- **magnetic_column** — 
- **temp_unit** — 
- **input_unit** — 
- **smooth_window** — As in ``curie_temperature_estimates``.
- **remove_holder** — As in ``curie_temperature_estimates``.
- **branches** — As in ``curie_temperature_estimates``.
- **method_kwargs** — As in ``curie_temperature_estimates``.
- **data_type** — As in ``curie_temperature_estimates``.
- **figsize** (`tuple`) — Figure size in inches (default (10, 10)).
- **return_figure** (`bool`) — If True, return ``(fig, axes)`` (default False).
- **save_path** (`str`) — If given, save the figure to this path.

**Returns**

- `(matplotlib.figure.Figure, numpy.ndarray of Axes) or None` — 

---

## `plot_day`

```python
plot_day(Mr, Ms, Bcr, Bc, Mr_Ms_lower=0.05, Mr_Ms_upper=0.5, Bc_Bcr_lower=1.5, Bc_Bcr_upper=4, plot_day_lines=True, plot_MD_slope=True, plot_SP_SD_mixing=[10, 15, 25, 30], plot_SD_MD_mixing=True, color='black', marker='o', label='sample', alpha=1, lc='black', lw=0.5, legend=True, figsize=(8, 6), show_plot=True, return_figure=True)
```

function to plot given Ms, Mr, Bc, Bcr values either as single values or list/array of values 
    plots Mr/Ms vs Bc/Bcr. 

**Parameters**

- **Ms** (`float or array - like`) — saturation magnetization
- **Mr** (`float or array - like`) — remanent magnetization
- **Bc** (`float or array - like`) — coercivity
- **Bcr** (`float or array - like`) — coercivity of remanence
- **color** (`str`) — color of the points. The default is 'black'.
- **marker** (`str`) — marker style of the points. The default is 'o'.
- **label** (`str`) — label for the points. The default is 'sample'.
- **alpha** (`float`) — transparency of the points. The default is 1.
- **lc** (`str`) — color of the lines. The default is 'black'.
- **lw** (`float`) — line width of the lines. The default is 0.5.
- **legend** (`bool`) — whether to show the legend. The default is True.
- **figsize** (`tuple`) — size of the figure. The default is (6,6).
- **show_plot** (`bool`) — whether to show the plot. The default is True.
- **return_figure** (`bool`) — whether to return the figure and axes objects. The default is True, so that a different function (plot_day_MagIC) can use it.

**Returns**

- `tuple or None` — - If return_figure is True (default), returns (fig, ax). - Otherwise, returns None.

---

## `plot_day_magic`

```python
plot_day_magic(specimen_data, by='specimen', Mr='hyst_mr_mass', Ms='hyst_ms_mass', Bcr='rem_bcr', Bc='hyst_bc', kwargs={})
```

Function to plot a Day plot from a MagIC specimens table.

**Parameters**

- **specimen_data** (`pandas.DataFrame`) — DataFrame containing the specimens data.
- **by** (`str`) — Column name to group by (default is 'specimen').
- **Mr** (`str`) — Column name for the remanence (default is 'hyst_mr_mass').
- **Ms** (`str`) — Column name for the saturation magnetization (default is 'hyst_ms_mass').
- **Bcr** (`str`) — Column name for the coercivity (default is 'hyst_bcr').
- **Bc** (`str`) — Column name for the coercivity of remanence (default is 'hyst_bc').
- ****kwargs** (`keyword arguments`) — Additional arguments to pass to the plotting function.

**Returns**

- `matplotlib.axes.Axes` — The axes object containing the plot.

---

## `plot_hyst_loop`

```python
plot_hyst_loop(field, magnetization, specimen_name, p=None, interactive=True, show_plot=True, return_figure=False, line_color='grey', line_width=1, label='', legend_location='bottom_right')
```

function to plot a hysteresis loop

**Parameters**

- **field** (`numpy array or list`) — hysteresis loop field values
- **magnetization** (`numpy array or list`) — hysteresis loop magnetization values

**Returns**

- `bokeh.plotting.figure` — 

---

## `plot_mpms_ac`

```python
plot_mpms_ac(experiment, frequency=None, phase='in', figsize=(6, 6), interactive=False, return_figure=False, show_plot=True, legend_location='upper left')
```

Plot AC susceptibility data from MPMS-X, optionally as interactive Bokeh.

**Parameters**

- **experiment** (`pandas.DataFrame`) — The experiment table from the MagIC contribution.
- **frequency** (`float or None`) — Frequency of AC measurement in Hz; None plots all frequencies.
- **phase** (`('in', out, both)`) — Which phase to plot.
- **figsize** (`tuple of float`) — Figure size for Matplotlib (width, height).
- **interactive** (`bool`) — If True, render with Bokeh for interactive exploration.
- **return_figure** (`bool`) — If True, return the figure object(s).
- **show_plot** (`bool`) — If True, display the plot.
- **legend_location** (`str`) — Location of the legend in Matplotlib terms.

**Returns**

- `fig, ax or (fig, axes) or Bokeh layout or None` — 

---

## `plot_mpms_dc`

```python
plot_mpms_dc(fc_data=None, zfc_data=None, rtsirm_cool_data=None, rtsirm_warm_data=None, fc_color='#1f77b4', zfc_color='#ff7f0e', rtsirm_cool_color='#17becf', rtsirm_warm_color='#d62728', fc_marker='d', zfc_marker='p', rtsirm_cool_marker='s', rtsirm_warm_marker='o', symbol_size=4, interactive=False, plot_derivative=False, return_figure=False, show_plot=True, drop_first=False, drop_last=False)
```

Plots MPMS DC data and optional derivatives, omitting empty panels.

Parameters:
    fc_data (DataFrame or None): Field-cooled data.
    zfc_data (DataFrame or None): Zero-field-cooled data.
    rtsirm_cool_data (DataFrame or None): RTSIRM cooling data.
    rtsirm_warm_data (DataFrame or None): RTSIRM warming data.
    fc_color, zfc_color, rtsirm_cool_color, rtsirm_warm_color (str):
        HEX color codes.
    fc_marker, zfc_marker, rtsirm_cool_marker, rtsirm_warm_marker (str):
        Matplotlib-style markers.
    symbol_size (int): Marker size.
    interactive (bool): If True, use Bokeh.
    plot_derivative (bool): If True, plot dM/dT curves.
    return_figure (bool): If True, return the figure/grid.
    show_plot (bool): If True, display the plot.
    drop_first (bool): If True, drop first row of each series.
    drop_last (bool): If True, drop last row of each series.

Returns:
    Bokeh grid or Matplotlib fig/axes tuple, or None.

---

## `plot_mpms_dc_interactive`

```python
plot_mpms_dc_interactive(measurements)
```

Create a UI to select specimen and plot MPMS data in Matplotlib or Bokeh.

Parameters:
    measurements (pandas.DataFrame): DataFrame with 'specimen' and
        'method_codes'.

---

## `plot_neel`

```python
plot_neel(Mr, Ms, Bc, color='black', marker='o', label='sample', alpha=1, lc='black', lw=0.5, legend=True, axis_scale='linear', figsize=(5, 5))
```

Generate a Néel plot (squareness-coercivity) of Mr/Ms versus Bc from hysteresis data.

This plot shows the ratio of remanent to saturation magnetization
(Mr/Ms) plotted against the coercivity (Bc). It is useful for 
characterizing magnetic domain states in rock magnetic samples.

**Parameters**

- **Mr** (`array - like`) — Saturation remanence values of the samples.
- **Ms** (`array - like`) — Saturation magnetization values of the samples.
- **Bc** (`array - like`) — Coercivity values of the samples.
- **color** (`str`) — Color of the scatter points. Default is "black".
- **marker** (`str`) — Marker style for scatter points. Default is "o".
- **label** (`str`) — Label for the sample to be displayed in the legend. Default is "sample".
- **alpha** (`float`) — Transparency of the scatter points. Default is 1 (opaque).
- **lc** (`str`) — Color of the grid lines. Default is "black".
- **lw** (`float`) — Line width of the grid lines. Default is 0.5.
- **legend** (`bool`) — Whether to show the legend. Default is True.
- **axis_scale** (`str`) — Scale for both axes: "linear" or "log". Default is "linear".
- **figsize** (`tuple of int`) — Figure size in inches (width, height). Default is (5, 5).

**Returns**

- `matplotlib.axes.Axes` — The matplotlib axes object containing the plot.

---

## `plot_neel_magic`

```python
plot_neel_magic(specimen_data, by='specimen', Mr='hyst_mr_mass', Ms='hyst_ms_mass', Bcr='rem_bcr', Bc='hyst_bc', kwargs={})
```

Function to plot a Day plot from a MagIC specimens table.

**Parameters**

- **specimen_data** (`pandas.DataFrame`) — DataFrame containing the specimens data.
- **by** (`str`) — Column name to group by (default is 'specimen').
- **Mr** (`str`) — Column name for the remanence (default is 'hyst_mr_mass').
- **Ms** (`str`) — Column name for the saturation magnetization (default is 'hyst_ms_mass').
- **Bcr** (`str`) — Column name for the coercivity (default is 'hyst_bcr').
- **Bc** (`str`) — Column name for the coercivity of remanence (default is 'hyst_bc').
- ****kwargs** (`keyword arguments`) — Additional arguments to pass to the plotting function.

**Returns**

- `matplotlib.axes.Axes` — The axes object containing the plot.

---

## `plot_unmixing_ensemble`

```python
plot_unmixing_ensemble(result, space='spectrum', n_draws=200, n_grid=300, show_components=True, figsize=(8, 5), title=None, colors=None, random_seed=None)
```

Overlay many draws of the model curves to visualize decomposition spread.

Rather than a single best fit with an error band, this draws the total
model and (optionally) the individual components for many bootstrap or
posterior samples, so the full range of decompositions consistent with
the data is visible directly -- including cases where components exchange
coercivity or amplitude between draws.

**Parameters**

- **result** (`dict`) — Result from unmixing_bootstrap or unmix_coercivity_bayes.
- **space** (`str`) — 'spectrum' (dM/dlog10 B) or 'curve' (measurement space).
- **n_draws** (`int`) — Number of draws to overlay (capped at the number available).
- **n_grid** (`int`) — Number of field points for the smooth curves.
- **show_components** (`bool`) — Overlay per-component curves in addition to the total.
- **figsize** (`tuple`) — Figure size.
- **title** (`str`) — Figure title.
- **colors** (`list`) — Per-component colors.
- **random_seed** (`None, int, or numpy.random.Generator`) — Seed for choosing which draws to plot.

**Returns**

- `tuple` — (fig, ax).

---

## `plot_unmixing_multistart`

```python
plot_unmixing_multistart(result, max_solutions=6, space='spectrum', n_grid=300, figsize=None, colors=None, marker_scale='uniform')
```

Visualize the distinct solutions found by a multi-start analysis.

Produces a panel of small multiples, one per distinct solution (ordered
best-fit first), each showing that solution's decomposition against the
data, plus a final parameter-space map that places every solution's
components on a coercivity-versus-proportion plot. Together these make
the non-uniqueness of the decomposition concrete: how many genuinely
different solutions the data admit and how each partitions the spectrum.

**Parameters**

- **result** (`dict`) — Result from unmixing_multistart (carries a 'multistart' entry).
- **max_solutions** (`int`) — Maximum number of distinct solutions to draw as small multiples (the best-fitting solutions are shown).
- **space** (`str`) — 'spectrum' (dM/dlog10 B) or 'curve' (measurement space) for the decomposition panels.
- **n_grid** (`int`) — Number of field points for the smooth model curves.
- **figsize** (`tuple`) — Figure size; a default is chosen from the panel count.
- **colors** (`list`) — Per-solution colors; defaults to the tab10 cycle.
- **marker_scale** (`str`) — How the solution-map markers are sized: 'uniform' (default, all equal, so no solution is visually privileged), 'n_hits' (size scaled into a bounded range by the number of starts that reached each solution), or 'weight' (size scaled by Akaike weight). On low-noise data the Akaike weight collapses onto the lowest-RSS solution, so 'uniform' or 'n_hits' better reflect that the solutions are alternatives; 'weight' is informative mainly when the noise is large enough to spread support across solutions. The number of starts that reached each solution is annotated on the map in every case.

**Returns**

- `tuple` — (fig, axes).

---

## `plot_unmixing_posterior`

```python
plot_unmixing_posterior(result, quantity='B_mean_mT', bins=40, figsize=None, colors=None)
```

Plot marginal uncertainty distributions of a component quantity.

Draws one histogram per component of the requested derived quantity from
the bootstrap or Bayesian draws, with the median and 95% interval marked.
This visualizes how tightly each component parameter is constrained,
including asymmetric and multimodal uncertainties that a single
standard-error value cannot convey.

**Parameters**

- **result** (`dict`) — Result from unmixing_bootstrap or unmix_coercivity_bayes.
- **quantity** (`str`) — Name of the quantity to plot (e.g. 'B_mean_mT', 'proportion', 'sd_log', 'location', 'dp', 'skew'). Must be present in the draws.
- **bins** (`int`) — Number of histogram bins.
- **figsize** (`tuple`) — Figure size; defaults to (7, 2.2 * n_components).
- **colors** (`list`) — Per-component colors; defaults to the matplotlib C0, C1, ... cycle.

**Returns**

- `tuple` — (fig, axes).

---

## `plot_unmixing_tradeoff`

```python
plot_unmixing_tradeoff(result, x='B_mean_mT', y='proportion', component=None, figsize=(5, 5), colors=None)
```

Scatter two component quantities across draws to reveal parameter trade-offs.

Overlapping coercivity components trade parameters against one another;
plotting one quantity against another across the bootstrap or posterior
draws exposes these correlations (and any multimodality) that marginal
intervals hide. By default every component is shown; pass a component
index to isolate one.

**Parameters**

- **result** (`dict`) — Result from unmixing_bootstrap or unmix_coercivity_bayes.
- **x** (`str`) — Quantity names for the two axes (e.g. 'B_mean_mT', 'proportion', 'sd_log').
- **y** (`str`) — Quantity names for the two axes (e.g. 'B_mean_mT', 'proportion', 'sd_log').
- **component** (`int`) — 1-based component index to plot alone; if None, all components are overlaid.
- **figsize** (`tuple`) — Figure size.
- **colors** (`list`) — Per-component colors.

**Returns**

- `tuple` — (fig, ax).

---

## `prepare_thermomag_branches`

```python
prepare_thermomag_branches(experiment, temperature_column='meas_temp', magnetic_column='susc_chi_mass', temp_unit='C', input_unit='K', smooth_window=0, remove_holder=True, window_type='hanning')
```

Preprocess a thermomagnetic experiment into clean heating/cooling branches.

This is the shared preprocessing step for the Curie temperature estimators
and thermomagnetic plots. It splits the measurement sequence into heating
and cooling branches (``split_heating_cooling``), converts temperatures to the
requested unit, optionally subtracts the per-branch minimum as an estimate
of the sample-holder background, sorts each branch by ascending
temperature, and optionally smooths each branch with an x-space moving
window (``smooth_moving_average``).

Subtracting the per-branch minimum assumes that the magnetic signal decays
to the holder background at the highest temperatures (i.e., the experiment
passes above the Curie temperature of all ferromagnetic phases). When that
assumption does not hold (e.g., a run that ends below the Curie
temperature), set ``remove_holder=False``.

**Parameters**

- **experiment** (`pandas.DataFrame`) — MagIC-formatted experiment DataFrame (rows in measurement order).
- **temperature_column** (`str`) — Name of the temperature column (default 'meas_temp').
- **magnetic_column** (`str`) — Name of the magnetization/susceptibility column (default 'susc_chi_mass').
- **temp_unit** (`(C, K)`) — Unit for the returned temperatures (default 'C').
- **input_unit** (`(K, C)`) — Unit of the temperatures in ``experiment`` (default 'K', the MagIC convention for ``meas_temp``).
- **smooth_window** (`float`) — Width of the smoothing window in units of ``temp_unit``. If 0 (default), no smoothing is applied and the smoothed arrays equal the raw arrays.
- **remove_holder** (`bool`) — Subtract the per-branch minimum value from each branch (default True).
- **window_type** (`(flat, hanning, hamming, bartlett, blackman)`) — Weighting function applied within each smoothing window by ``smooth_moving_average`` (default 'hanning'). Only used when ``smooth_window > 0``. The options are:  - 'flat': uniform weights, i.e. a simple unweighted running mean. - 'hanning': raised-cosine (Hann) taper; weights fall smoothly to   zero at the window edges. A good general-purpose default that   suppresses edge/ringing artifacts. - 'hamming': raised-cosine taper similar to 'hanning' but with   nonzero end weights, giving slightly less edge attenuation. - 'bartlett': triangular taper; weights decrease linearly from the   window center to zero at the edges. - 'blackman': three-term cosine taper that is more strongly peaked   than 'hanning'/'hamming', giving the heaviest smoothing (widest   effective averaging) of the tapered options.  All options other than 'flat' are the correspondingly named ``numpy`` window functions.

**Returns**

- `dict` — ``{'heating': branch or None, 'cooling': branch or None}`` where each branch is a dict with keys ``'T'`` and ``'y'`` (smoothed arrays, ascending temperature) and ``'raw_T'`` and ``'raw_y'`` (unsmoothed arrays, ascending temperature). A branch with no measurements is ``None``.

---

## `process_backfield_data`

```python
process_backfield_data(experiment, field='treat_dc_field', magnetization='magn_mass', smooth_mode='lowess', smooth_frac=0.0, drop_first=False)
```

Function to process the backfield data including shifting the magnetic 
moment to be positive values taking the log base 10 of the magnetic 
field values and writing these new fields into the experiment attribute 
table

**Parameters**

- **experiment** (`DataFrame`) — DataFrame containing the backfield data
- **field** (`str`) — The name of the treatment field column in the DataFrame
- **magnetization** (`str`) — The name of the magnetization column in the DataFrame
- **smooth_mode** (`str`) — The smoothing mode to be used, either 'lowess' or 'spline'
- **smooth_frac** (`float`) — Fraction of the data to be used for LOWESS smoothing, value must be between 0 and 1
- **drop_first** (`bool`) — Whether to drop the first data point or not in some cases you may want to drop the first data point to avoid negative log values

**Returns**

- `DataFrame` — The processed experiment DataFrame with new attributes.

---

## `process_hyst_loop`

```python
process_hyst_loop(field, magnetization, specimen_name='', show_results_table=True, show_plot=True, NL_fit=False, centering_protocol='legacy', fit_open_loop=False, fit_linear_loop=False)
```

Process a magnetic hysteresis loop using the IRM decision tree workflow.

This function performs a complete analysis of a hysteresis loop, including gridding, centering, drift correction,
high-field correction, and extraction of key magnetic parameters. The workflow follows best practices in rock magnetism
and outputs both a summary of results and a Bokeh plot visualizing the various processing steps.

The inputs need not come from a MagIC measurements table: any pair of
field and magnetization sequences (lists, arrays, or dataframe columns)
can be processed. Inputs are passed through `sanitize_hyst_inputs`,
so non-finite measurement pairs are dropped with a report, numeric
strings are converted, and either field sweep order (starting from
positive or negative saturation) is accepted.

Two decision-tree exits terminate processing early, in both cases
returning the full result key set with the undefined quantities reported
as NaN/None so batch tables keep a stable schema:

- a loop that is statistically linear (whole-loop lack-of-fit test) is
  dominated by paramagnetic or diamagnetic material; only the high-field
  susceptibility (from the whole-loop regression) is reported. Passing
  fit_linear_loop=True overrides this exit and processes the loop in
  full;
- a loop that remains open at the highest fields (closure test) contains
  unsaturated high-coercivity phases, so Ms and chi_HF cannot be
  separated; the slope-independent parameters (Mr and Brh, computed from
  Mrh in which linear-in-field contributions cancel) and the data
  quality statistics are reported. Passing fit_open_loop=True overrides
  this exit and proceeds with the high-field fitting.

**Parameters**

- **field** (`array_like`) — Array of applied magnetic field values in tesla (the chi_HF unit conversions assume tesla; a warning is printed if the values appear to be in mT or Oe).
- **magnetization** (`array_like`) — Array of magnetization values (same length as `field`), in any consistent unit; mass-normalized Am2/kg matches MagIC conventions.
- **specimen_name** (`str`) — Identifier for the specimen, used for labeling plots.
- **show_results_table** (`bool`) — If True (default), display a summary table of key parameters using Bokeh.
- **show_plot** (`bool`) — If True (default), display the Bokeh plot of the hysteresis loop and processing steps.
- **NL_fit** (`bool`) — If True, force non-linear high-field fitting regardless of the saturation test result (default is False). Because the approach-to-saturation fit exists precisely for unsaturated loops, NL_fit=True also proceeds through the open-loop exit (it implies fit_open_loop=True).
- **centering_protocol** (`(legacy, iterative)`) — Centering workflow to apply before drift and high-field corrections. Defaults to 'legacy' for backward compatibility.
- **fit_open_loop** (`bool`) — If True, proceed with the high-field fitting (and the Ms estimate) even when the closure test flags the loop as open. Default False: open loops exit with the slope-independent parameters and data quality statistics, since Ms and chi_HF cannot be separated for an unsaturated loop. NL_fit=True implies this behavior. Note that residual instrument drift can leave a spurious positive high-field Mrh signal that trips the closure test on a visually closed loop (particularly for loops measured from negative saturation); inspect the loop and pass fit_open_loop=True in such cases.
- **fit_linear_loop** (`bool`) — If True, process a statistically linear loop in full rather than terminating with chi_HF only (default False). Useful when a weak ferromagnetic signal near the noise level is of interest despite the loop passing the whole-loop linearity test; the ferromagnetic parameters from such a loop should be interpreted alongside the quality statistics.

**Returns**

- `dict` — Dictionary containing the following keys:     - 'gridded_H': gridded field values     - 'gridded_M': gridded magnetization values     - 'linearity_test_results': results of the initial linearity test     - 'loop_is_linear': whether the loop passes the linearity test     - 'FNL': F statistic for whole-loop nonlinearity (lack-of-fit F ratio)     - 'loop_centering_results': results of centering optimization     - 'centered_H': centered field values     - 'centered_M': centered magnetization values     - 'drift_corrected_M': drift-corrected magnetization     - 'slope_corrected_M': slope-corrected magnetization     - 'loop_closure_test_results': results of closure test     - 'loop_is_closed': whether the loop is closed     - 'loop_saturation_stats': saturation test results     - 'loop_is_saturated': whether the loop is saturated     - 'M_sn', 'Q': quality metrics from centering     - 'H', 'Mr', 'Mrh', 'Mih', 'Me', 'Brh': characteristic field and moment parameters     - 'sigma': shape parameter (Fabian, 2003)     - 'chi_HF': high-field susceptibility     - 'FNL60', 'FNL70', 'FNL80': high-field nonlinearity F statistics for windows       starting at 60%, 70%, and 80% of the maximum field     - 'Ms': saturation magnetization     - 'Bc': coercive field     - 'M_sn_f', 'Qf': quality metrics for ferromagnetic component     - 'Fnl_lin': F statistic for improvement of the nonlinear over the linear       high-field fit (None if the loop is saturated and no nonlinear fit is made)     - 'plot': Bokeh figure with overlaid processing steps

---

## `process_hyst_loops`

```python
process_hyst_loops(hyst_experiments, measurements, field_col='meas_field_dc', magn_col='magn_mass', show_results_table=True, show_plots=True, centering_protocol='legacy', fit_open_loop=False, fit_linear_loop=False)
```

Process multiple hysteresis loops in batch.

**Parameters**

- **hyst_experiments** (`DataFrame`) — Must contain columns "experiment" and "specimen".
- **measurements** (`DataFrame`) — Must contain an "experiment" column and the data columns.
- **field_col** (`str`) — Name of the column in `measurements` holding field values. Defaults to "meas_field_dc".
- **magn_col** (`str`) — Name of the column in `measurements` holding magnetization values. Defaults to "magn_mass".
- **show_results_table** (`bool`) — If True, display the summary table below each plot.
- **show_plots** (`bool`) — If True, display the hysteresis plots for each specimen.
- **centering_protocol** (`(legacy, iterative)`) — Centering workflow to pass through to process_hyst_loop. Defaults to 'legacy' for backward compatibility.
- **fit_open_loop** (`bool`) — Passed through to process_hyst_loop: if True, high-field fitting proceeds even for loops the closure test flags as open (default False).
- **fit_linear_loop** (`bool`) — Passed through to process_hyst_loop: if True, statistically linear loops are processed in full rather than terminating with chi_HF only (default False).

**Returns**

- `pandas.DataFrame` — DataFrame with hysteresis results for each experiment. Has a numeric index with 'specimen' and 'experiment' as columns.

---

## `prorated_drift_correction`

```python
prorated_drift_correction(field, magnetization, descending_first=True)
```

function to correct for the linear drift of a hysteresis loop
    take the difference between the magnetization measured at the maximum field on the upper and lower branches
    apply linearly prorated correction of M(H)
    this should be applied to the gridded data

The prorated ramp runs in measurement-time order, and the arrays are
expected in canonical order (descending upper branch first, as produced
by `grid_hyst_loop`). For a loop originally measured from negative
saturation, pass descending_first=False (detected from the raw field
values with `measured_descending_first`) so the ramp is applied in true
time order rather than with the opposite time sense.

**Parameters**

- **field** (`numpy array`) — field values, in canonical (descending-upper-branch-first) order
- **magnetization** (`numpy array`) — magnetization values
- **descending_first** (`bool`) — Whether the loop was originally measured with the descending branch first (default True). Use `measured_descending_first` on the raw field values to determine this for a gridded loop.

**Returns**

- `numpy array` — corrected magnetization values

---

## `register_unmixing_method`

```python
register_unmixing_method(name, function)
```

Register a custom coercivity unmixing method.

Registered methods become available through unmix_coercivity and
unmix_backfield_experiments alongside the built-in 'spectrum', 'curve',
and 'maxunmix' methods. The function must accept
(x, magnetization, n_components=None, initial_parameters=None,
curve_type='backfield', vary_skew=True, **kwargs) and return the
standardized result dictionary (see unmix_coercivity_spectrum); the
component model helpers (skewnormal_pdf, coercivity_spectrum_model,
coercivity_curve_model, ...) can be reused when implementing new
methods.

**Parameters**

- **name** (`str`) — Name under which the method is registered.
- **function** (`callable`) — The method implementation.

**Returns**

- `None` — 

---

## `sanitize_hyst_inputs`

```python
sanitize_hyst_inputs(field, magnetization, drop_nonfinite=True)
```

Coerce hysteresis loop inputs into clean float arrays for processing.

This helper makes the processing functions usable with data from any
source: MagIC measurement table columns (including ones read as text),
plain Python lists, or arrays exported from instrument software.

**Parameters**

- **field** (`array_like`) — Applied field values. Expected in tesla: the chi_HF unit conversions in `linear_HF_fit`, `hyst_slope_correction`, and the nonlinear fits assume tesla, so a warning is printed if the values appear to be in mT or Oe (max |field| > 20).
- **magnetization** (`array_like`) — Magnetization or moment values, in any unit that is consistent across the loop (mass-normalized Am2/kg matches MagIC conventions).
- **drop_nonfinite** (`bool`) — If True (default), measurement pairs where either value is NaN or infinite are dropped with a printed report. If False, a ValueError is raised when non-finite values are present.

**Returns**

- `field, magnetization : numpy.ndarray` — Equal-length float arrays with only finite values.

---

## `select_n_components`

```python
select_n_components(x, magnetization, method=DEFAULT_UNMIX_METHOD, min_components=1, max_components=4, criterion='parsimony', min_improvement=0.02, noise_level=None, reduced_chi2_target=1.0, curve_type='backfield', vary_skew=DEFAULT_UNMIX_VARY_SKEW, verbose=False, kwargs={})
```

Choose the number of coercivity components by a parsimony rule.

Fits models with a range of component counts and selects the simplest
one that adequately describes the data. This is deliberately different
from minimizing an information criterion or maximizing the Bayesian
evidence: with high-resolution, low-noise curves those measures tend to
keep favoring more components indefinitely, because a real coercivity
distribution is never exactly log-Gaussian and each added component
removes a little more systematic misfit. Following the parsimony
principle emphasized by Egli (2003) and Heslop (2015), an extra
component is accepted only when it produces a large enough improvement,
so a simpler adequate model is preferred.

Two selection criteria are provided:

- 'parsimony' (default): an added component is retained only if it
  reduces the residual sum of squares by at least `min_improvement`
  times the baseline (single-component) residual. Because the second
  component typically removes most of the baseline misfit while a
  spurious third component removes only a tiny fraction of it, this
  robustly stops at the mineralogically meaningful count regardless of
  the noise level.
- 'chi2': the simplest model whose reduced chi-square (using
  `noise_level`, estimated with estimate_measurement_noise if not
  given) falls at or below `reduced_chi2_target`, i.e. the simplest
  model that fits to within the measurement noise. The noise estimator
  is spacing-aware and unbiased on the coarse, log-spaced field grids of
  backfield/IRM data, but 'chi2' remains the less robust criterion: like
  the information criteria and the Bayesian evidence it tends to
  over-select on high-resolution, low-noise curves (once the noise is
  not over-estimated, any small departure of the data from a log-Gaussian
  pushes the reduced chi-square of a real-count model just above one), and
  it is sensitive to the residual scatter of the noise estimate. Prefer
  'parsimony', or supply a trusted `noise_level` and a `reduced_chi2_target`
  slightly above 1, when using it. For spectrum-space methods ('spectrum',
  'maxunmix', or a Bayesian fit with space='spectrum') the residuals and
  hence the noise are in dM/dlog10(B) units, so the noise is estimated on
  the finite-difference spectrum rather than the measured curve; that
  spectrum noise is mildly correlated, which the estimator does not model,
  making the spectrum-space chi2 more approximate still.

In all cases the returned table reports the fit statistics, the
fractional RSS improvement per added component, and (when a noise level
is available) the reduced chi-square, so the selection can be inspected
and overridden.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **magnetization** (`array - like`) — Remanence curve values at x.
- **method** (`str`) — Registered unmixing method used for each fit (default DEFAULT_UNMIX_METHOD, i.e. 'spectrum').
- **min_components** (`int`) — Range of component counts to consider.
- **max_components** (`int`) — Range of component counts to consider.
- **criterion** (`str`) — 'parsimony' or 'chi2'.
- **min_improvement** (`float`) — For 'parsimony': minimum fraction of the baseline residual an added component must explain to be retained (default 0.02).
- **noise_level** (`float`) — Measurement noise standard deviation for 'chi2'; estimated from the data if not given.
- **reduced_chi2_target** (`float`) — For 'chi2': the reduced chi-square at or below which a model is considered adequate (default 1.0).
- **curve_type** (`str`) — 'backfield' or 'acquisition'.
- **vary_skew** (`bool`) — Whether skew varies during fitting.
- **verbose** (`bool`) — Print the selection outcome.
- ****kwargs** — Passed to the unmixing method.

**Returns**

- `tuple` — (selected_n, table, results) where selected_n is the chosen number of components, table is a DataFrame of per-model statistics (with a boolean 'selected' column), and results maps component count -> result dictionary.

---

## `skewnormal_cdf`

```python
skewnormal_cdf(x, location, dp, skew=0.0)
```

Skew-normal cumulative distribution function.

Evaluated analytically as Phi(z) - 2*T(z, skew) where T is Owen's T
function (scipy.special.owens_t).

**Parameters**

- **x** (`array - like`) — Points at which to evaluate the CDF (log10 of field in mT).
- **location** (`float`) — Location parameter (log10 mT).
- **dp** (`float`) — Scale parameter (log10 units); must be positive.
- **skew** (`float`) — Shape parameter alpha (default 0).

**Returns**

- `numpy.ndarray` — CDF values between 0 and 1.

---

## `skewnormal_pdf`

```python
skewnormal_pdf(x, location, dp, skew=0.0)
```

Skew-normal probability density function (unit area).

With skew = 0 this is a Gaussian with mean = location and standard
deviation = dp; in log10(B) coordinates that Gaussian is the log-Gaussian
coercivity distribution of Robertson & France (1994). Nonzero skew
follows the Azzalini (1985) formulation: negative values skew the
distribution toward low values (tail to the left).

**Parameters**

- **x** (`array - like`) — Points at which to evaluate the density (log10 of field in mT).
- **location** (`float`) — Location parameter (log10 mT). Equal to the mean only when skew = 0.
- **dp** (`float`) — Scale parameter (log10 units); must be positive. Equal to the standard deviation only when skew = 0.
- **skew** (`float`) — Shape parameter alpha of the Azzalini skew-normal (default 0).

**Returns**

- `numpy.ndarray` — Density values with unit integrated area.

---

## `skewnormal_stats`

```python
skewnormal_stats(location, dp, skew=0.0)
```

Moments and characteristic points of a skew-normal distribution.

**Parameters**

- **location** (`float`) — Location parameter (log10 mT).
- **dp** (`float`) — Scale parameter (log10 units).
- **skew** (`float`) — Shape parameter alpha.

**Returns**

- `dict` — With keys 'mean', 'std', 'median', and 'mode', all in the same (log10) units as location and dp. For skew = 0 all of mean, median, and mode equal location and std equals dp.

---

## `smooth_moving_average`

```python
smooth_moving_average(x, y, x_window, window_type='hanning', pad_mode='edge', return_variance=False)
```

Smooth y vs x using an x-space moving window and numpy window functions.

Parameters:
    x (array-like):
        1-D sequence of independent variable values.
    y (array-like):
        1-D sequence of dependent variable values.
    x_window (float):
        Width of the x-window centered on each point; must be >= 0.
        If zero, no smoothing is applied.
    window_type (str, optional):
        One of ['flat', 'hanning', 'hamming', 'bartlett', 'blackman'].
        'flat' is a simple running mean. Defaults to 'hanning'.
    pad_mode (str, optional):
        Mode for numpy.pad to reduce edge artifacts (e.g., 'edge',
        'constant', 'nearest'). Defaults to 'edge'.
    return_variance (bool, optional):
        If True, return weighted variances of x and y as well.
        Otherwise, only return smoothed x and y. Defaults to False.

Returns:
    tuple: ``(smoothed_x, smoothed_y)`` by default, or
    ``(smoothed_x, smoothed_y, x_var, y_var)`` when return_variance is
    True. ``smoothed_x`` and ``smoothed_y`` are the window-averaged arrays
    (same length as the inputs); ``x_var`` and ``y_var`` are the
    corresponding per-point weighted variances within each window (in the
    squared units of x and y), a measure of local spread. When
    x_window is 0 the inputs are returned unchanged and the variances are
    zero.

---

## `specimen_experiment_selection_interactive`

```python
specimen_experiment_selection_interactive(measurements)
```

Creates interactive dropdown widgets for selecting a specimen and its associated
experiment from a measurements DataFrame.

**Parameters**

- **measurements** (`pd.DataFrame`) — DataFrame containing measurement data with at least two columns: 'specimen' and 'experiment'. The 'specimen' column holds the specimen names while the 'experiment' column holds the experiment identifiers associated with each specimen.

**Returns**

- `tuple of ipywidgets.Dropdown` — A tuple containing two dropdown widgets. The first widget allows for selecting a specimen, and the second widget allows for selecting an experiment associated with the chosen specimen. The experiment dropdown is dynamically updated based on the specimen selection.

---

## `specimen_selection_interactive`

```python
specimen_selection_interactive(measurements)
```

Creates and displays a dropdown widget for selecting a specimen from a given
DataFrame of measurements.

**Parameters**

- **measurements** (`pd.DataFrame`) — The DataFrame containing measurement data with a column 'specimen'. It is expected to have at least this column where 'specimen' identifies the specimen name.

**Returns**

- `ipywidgets.Dropdown` — A dropdown widget allowing for the selection of a specimen. The initial selection in the dropdown is set to the first specimen option.

---

## `split_heating_cooling`

```python
split_heating_cooling(experiment, temperature_column='meas_temp', magnetic_column='susc_chi_mass')
```

Split a thermomagnetic curve into heating and cooling portions.

The sequence is split at the temperature turning point (the global
maximum): measurements up to and including the peak form the heating
branch and the descending remainder forms the cooling branch. A run whose
temperature never descends after its peak returns an empty cooling branch,
and a run that descends from its first measurement returns an empty heating
branch. Rows with non-finite temperature or magnetic values are dropped.

Splitting at the turning point rather than on the sign of each local step
keeps repeated furnace-stabilization temperatures (where the step is zero)
and noisy readings on the heating ramp within the heating branch. A
point-by-point classification instead misroutes those points into a
spurious cooling branch, so a heating-only run with duplicated or noisy
temperatures would otherwise produce a phantom cooling curve. This assumes
a single heat-then-cool trajectory, the standard thermomagnetic protocol.

**Parameters**

- **experiment** (`pandas.DataFrame`) — the experiment data (rows in measurement order)
- **temperature_column** (`str`) — name of the temperature column (default 'meas_temp')
- **magnetic_column** (`str`) — name of the magnetization/susceptibility column (default 'susc_chi_mass')

**Returns**

- `numpy.ndarray` — temperatures for the heating cycle (measurement order)
- `numpy.ndarray` — magnetization/susceptibility for the heating cycle
- `numpy.ndarray` — temperatures for the cooling cycle (measurement order)
- `numpy.ndarray` — magnetization/susceptibility for the cooling cycle

---

## `split_hyst_loop`

```python
split_hyst_loop(field, magnetization)
```

function to split a hysteresis loop into upper and lower branches
    at the reversal of the applied field sweep

The loop reversal is located with `find_hyst_turning_point`, which
tolerates repeated field plateaus and minor field glitches. Loops measured
in either sweep order are supported: the branch measured from the positive
field extreme downward is returned as the upper branch regardless of
whether it was measured first or second.

**Parameters**

- **field** (`numpy array or list`) — hysteresis loop field values
- **magnetization** (`numpy array or list`) — hysteresis loop magnetization values

**Returns**

- `list` — [field, magnetization] for the upper branch, in ascending field order
- `list` — [field, magnetization] for the lower branch, in ascending field order

---

## `symmetric_averaging_drift_correction`

```python
symmetric_averaging_drift_correction(field, magnetization)
```

Apply symmetric averaging drift correction to a hysteresis loop.

This function corrects drift in magnetic hysteresis loop data by averaging the upper branch and
the inverted lower branch of the magnetization curve, then adjusting for tip-to-tip separation.
The corrected magnetization is constructed by concatenating the reversed, drift-corrected upper branch
and its inverted counterpart, restoring symmetry to the loop.

**Parameters**

- **field** (`array_like`) — Array of applied magnetic field values for the hysteresis loop.
- **magnetization** (`array_like`) — Array of measured magnetization values corresponding to `field`.

**Returns**

- `numpy.ndarray` — Array of drift-corrected magnetization values, symmetrically constructed for the full loop. 

**Examples**

```python
>>> field = np.linspace(-1, 1, 200)
>>> magnetization = some_hysteresis_measurement(field)
>>> corrected = symmetric_averaging_drift_correction(field, magnetization)
```

---

## `thermomag_derivative`

```python
thermomag_derivative(temps, mags, drop_first=False, drop_last=False)
```

Calculates the derivative of magnetization with respect to temperature and optionally
drops the data corresponding to the highest and/or lowest temperature.

Parameters:
    temps (pd.Series): A pandas Series representing the temperatures at which
                       magnetization measurements were taken.
    mags (pd.Series): A pandas Series representing the magnetization measurements.
    drop_last (bool): Optional; whether to drop the last row from the resulting
                      DataFrame. Defaults to False. Useful when there is an
                      artifact associated with the end of the experiment.
    drop_first (bool): Optional; whether to drop the first row from the resulting
                       DataFrame. Defaults to False. Useful when there is an
                      artifact associated with the start of the experiment.

Returns:
    pd.DataFrame: A pandas DataFrame with two columns:
                  'T' - Midpoint temperatures for each temperature interval.
                  'dM_dT' - The derivative of magnetization with respect to temperature.
                  If drop_last is True, the last temperature point is excluded.
                  If drop_first is True, the first temperature point is excluded.

---

## `unmix_backfield_curve`

```python
unmix_backfield_curve(x, magnetization, n_components=None, initial_parameters=None, curve_type='backfield', vary_skew=DEFAULT_UNMIX_VARY_SKEW, fit_offset=True, weights=None, dp_bounds=(0.01, 2.0), skew_bounds=(-10.0, 10.0))
```

Unmix a remanence curve by fitting cumulative components directly.

Fits the measured curve M(log10 B) with a sum of skew-normal CDF
components (plus an optional constant offset), avoiding numerical
differentiation and smoothing entirely. The component parameterization
is identical to unmix_coercivity_spectrum, so results from the two
data spaces are directly comparable. Fitting in measurement space uses
the raw measurements with their original noise structure; fitting in
spectrum space can be more visually intuitive. Agreement between the
two approaches is a good indication of a robust unmixing model.

For backfield data processed with process_backfield_data, pass
x = 'log_dc_field' and magnetization = 'magn_mass_shift'. Note that the
shifted backfield curve spans twice the saturation remanence, so each
fitted 'contribution' is twice the remanence carried by that component;
'proportion' values are unaffected.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **magnetization** (`array - like`) — Remanence curve values at x (shifted positive for backfield data).
- **n_components** (`int`) — Number of components. Required if initial_parameters is None.
- **initial_parameters** (`pandas.DataFrame`) — Initial guesses (see unmix_coercivity_spectrum). If None, automatic estimates are derived from the finite-difference spectrum.
- **curve_type** (`str`) — 'backfield' (decaying curve) or 'acquisition' (growing curve).
- **vary_skew** (`bool`) — If False, skew values are fixed at their initial values.
- **fit_offset** (`bool`) — Whether to fit a constant baseline offset (default True).
- **weights** (`array - like`) — Multiplicative weights applied to the residuals.
- **dp_bounds** (`tuple`) — Bounds as in unmix_coercivity_spectrum.
- **skew_bounds** (`tuple`) — Bounds as in unmix_coercivity_spectrum.

**Returns**

- `dict` — Standardized result dictionary (see unmix_coercivity_spectrum), additionally including the fitted 'offset' and 'se_offset'.

---

## `unmix_backfield_experiments`

```python
unmix_backfield_experiments(measurements, experiments=None, n_components=2, method=DEFAULT_UNMIX_METHOD, initial_parameters=None, vary_skew=DEFAULT_UNMIX_VARY_SKEW, n_boot=0, resample='cases', proportion=1.0, noise_level=None, smooth_mode='spline', smooth_frac=0.0, drop_first=False, field='treat_dc_field', magnetization='magn_mass', random_seed=None, verbose=True, method_kwargs={})
```

Batch coercivity unmixing of backfield experiments in a MagIC
measurements table.

Each experiment is processed with process_backfield_data and unmixed
with the requested method (any name registered in UNMIXING_METHODS,
dispatched through unmix_coercivity). When initial parameters are not
supplied they are estimated automatically per experiment; supplying a
common initial-parameter table (or a per-experiment dict, e.g. built
with coercivity_unmixing_interactive) enforces a consistent starting
model across specimens, which aids comparability of the resulting
components.

**Parameters**

- **measurements** (`pandas.DataFrame`) — MagIC measurements table (must include 'experiment', 'specimen', 'method_codes', and the field/magnetization columns).
- **experiments** (`list`) — Experiment names to process. Defaults to all experiments whose method_codes include 'LP-BCR-BF'.
- **n_components** (`int`) — Number of components fit to each experiment (default 2).
- **method** (`str`) — Registered unmixing method name: 'spectrum', 'curve', 'maxunmix', or a custom method added with register_unmixing_method.
- **initial_parameters** (`pandas.DataFrame or dict`) — Either a single initial-parameter table applied to every experiment, or a dict mapping experiment name -> table. Experiments missing from the dict fall back to automatic estimation.
- **vary_skew** (`bool`) — Whether skew parameters vary during fitting.
- **n_boot** (`int`) — If > 0, ensure each result carries a bootstrap with this many replicates (methods that bootstrap internally, like 'maxunmix', are not re-bootstrapped).
- **resample** (`see unmixing_bootstrap.`) — 
- **proportion** (`see unmixing_bootstrap.`) — 
- **noise_level** (`see unmixing_bootstrap.`) — 
- **smooth_mode** (`see process_backfield_data.`) — The defaults (spline with smooth_frac=0) leave the data unsmoothed.
- **smooth_frac** (`see process_backfield_data.`) — The defaults (spline with smooth_frac=0) leave the data unsmoothed.
- **drop_first** (`see process_backfield_data.`) — The defaults (spline with smooth_frac=0) leave the data unsmoothed.
- **field** (`str`) — Column names in the measurements table.
- **magnetization** (`str`) — Column names in the measurements table.
- **random_seed** (`None, int, or numpy.random.Generator`) — Seed for reproducible bootstraps.
- **verbose** (`bool`) — Print progress and failures.
- ****method_kwargs** — Additional keyword arguments passed to the unmixing method.

**Returns**

- `tuple` — (components_df, results) where components_df is a tidy DataFrame with one row per experiment and component (parameters, derived coercivities, uncertainties, fit statistics, and Bcr) and results is a dict mapping experiment name -> full result dictionary.

---

## `unmix_coercivity`

```python
unmix_coercivity(x, magnetization, method=DEFAULT_UNMIX_METHOD, n_components=None, initial_parameters=None, curve_type='backfield', vary_skew=DEFAULT_UNMIX_VARY_SKEW, kwargs={})
```

Unmix a remanence curve into coercivity components with a named method.

This is the common entry point to the unmixing approaches implemented
in rockmagpy (and to any user-registered methods; see
register_unmixing_method). All methods take the measured remanence
curve -- not a precomputed derivative -- and share the same skew-normal
component parameterization, so their results are directly comparable.

Built-in methods:

- 'spectrum': fits the finite-difference coercivity spectrum
  dM/dlog10(B) with skew-normal components (Kruiver et al., 2001;
  Egli, 2003 lineage). Point estimates only; combine with
  unmixing_bootstrap for uncertainties.
- 'curve': fits the measured curve directly with cumulative
  (CDF) components, avoiding numerical differentiation entirely.
- 'maxunmix': the 'spectrum' fit plus the MAX UnMix resampling
  uncertainty scheme (Maxbauer et al., 2016): 95% case resampling with
  2% noise, 100 replicates by default.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT), e.g. 'log_dc_field' from process_backfield_data.
- **magnetization** (`array - like`) — Remanence curve values at x (e.g. 'magn_mass_shift').
- **method** (`str`) — Name of a registered unmixing method (default 'spectrum').
- **n_components** (`int`) — Number of components (required if initial_parameters is None).
- **initial_parameters** (`pandas.DataFrame`) — Initial guesses with columns 'contribution', 'location', 'dp', 'skew'; automatic estimates are used when omitted.
- **curve_type** (`str`) — 'backfield' or 'acquisition'.
- **vary_skew** (`bool`) — Whether skew parameters vary during fitting.
- ****kwargs** — Passed through to the method implementation (e.g. n_boot, proportion, param_noise for 'maxunmix'; fit_offset for 'curve'; dp_bounds, skew_bounds for the spectrum-based methods).

**Returns**

- `dict` — Standardized result dictionary (see unmix_coercivity_spectrum).

---

## `unmix_coercivity_bayes`

```python
unmix_coercivity_bayes(x, magnetization, n_components=2, curve_type='backfield', space='curve', vary_skew=False, fit_offset=True, priors=None, nlive=250, dlogz=0.1, sample='rslice', random_seed=None, n_grid=200, n_posterior_curves=300, verbose=False)
```

Bayesian coercivity unmixing by nested sampling (requires dynesty).

The remanence data are modeled as a sum of skew-normal components, with
the noise standard deviation treated as a free parameter, and sampled
with static nested sampling (Skilling, 2006) as implemented in dynesty
(Speagle, 2020), which also returns the Bayesian evidence (logz) -- the
principled criterion for choosing the number of components (compare logz
between runs with different n_components). The fit can be performed in
either of two data spaces (see the `space` argument):

- space='curve' (default): the measured curve M(B) is fit directly with
  cumulative skew-normal components plus a constant offset. An
  independent Gaussian noise model is defensible here, and no numerical
  differentiation is required. This is the more conservative choice.
- space='spectrum': the finite-difference coercivity spectrum
  dM/dlog10(B) is fit with skew-normal densities (no offset). Fitting the
  derivative directly reproduces the coercivity-distribution peak that a
  curve fit can under-represent, and lets skewness be constrained by the
  peak shape. The trade-off is that differencing correlates adjacent
  points, so the i.i.d. Gaussian likelihood used here is an approximation
  (the same one the least-squares and MAX UnMix spectrum fits make);
  credible intervals in this space should be read with that caveat.

The component parameterization (contribution=area, location, dp, skew) is
identical in both spaces, so their results are directly comparable, and
comparing them is a useful robustness check.

Unlike bootstrap resampling of a single fit, the posterior represents
the full range of component decompositions consistent with the data
and priors: parameter trade-offs between overlapping components appear
as wide, correlated, and possibly multimodal posterior distributions
rather than being hidden by a single optimizer solution.

Component locations are sampled as ordered order-statistics, which
fixes component labels without distorting the prior. Default priors
are weakly informative (locations uniform across the measured field
range, dispersions log-uniform on [0.02, 1.0] decades, contributions
uniform up to ~3x the data range); mineralogical knowledge can be
injected through the priors argument.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **magnetization** (`array - like`) — Remanence curve values at x (e.g. 'magn_mass_shift').
- **n_components** (`int`) — Number of components (default 2).
- **curve_type** (`str`) — 'backfield' or 'acquisition'.
- **space** (`str`) — 'curve' (fit the measured curve, default) or 'spectrum' (fit the finite-difference dM/dlog10(B) spectrum). See the summary above for the trade-offs; with space='spectrum' the offset is not used.
- **vary_skew** (`bool`) — Sample component skew (uniform prior on [-10, 10]); default False (symmetric log-Gaussian components). This deliberately deviates from DEFAULT_UNMIX_VARY_SKEW: free skew multiplies the nested-sampling cost and is usually better constrained through explicit `priors` windows (e.g. from mineral_priors) than left fully free.
- **fit_offset** (`bool`) — Include a constant baseline offset (default True; ignored when space='spectrum').
- **priors** (`dict`) — Overrides for the default prior bounds. Recognized keys: 'mean', 'location', 'dp', 'contribution', 'skew' map to a list of (low, high) tuples, one per component (in log10 units for 'mean'/'location'/'dp', magnetization units for 'contribution'); 'offset' and 'noise' map to a single (low, high) tuple in magnetization units. A 'mean' window constrains each component's MEAN coercivity (log10 mT): the mean is sampled uniformly in the window and the skew-normal location is derived from the sampled dp and skew, so the window means what it says even for skewed components (this is what mineral_priors produces). A 'location' window instead constrains the raw location parameter directly. Either replaces the weakly-informative ordered-uniform default, so the windows should be non-overlapping or ordered to keep component labels meaningful. 'mean' takes precedence over 'location' if both are given.
- **nlive** (`int`) — Number of live points (default 250).
- **dlogz** (`float`) — Evidence convergence tolerance (default 0.1).
- **sample** (`str`) — dynesty sampling method (default 'rslice'). Slice sampling is robust to the thin, curved likelihood ridges that overlapping components produce; the dynesty default uniform-ellipsoid sampler can stall on such geometries.
- **random_seed** (`None, int, or numpy.random.Generator`) — Seed for reproducibility.
- **n_grid** (`int`) — Grid size for posterior model bands.
- **n_posterior_curves** (`int`) — Number of posterior draws used for the model bands (default 300).
- **verbose** (`bool`) — Show dynesty progress.

**Returns**

- `dict` — Standardized result dictionary (parameters set to posterior medians) with an added 'bayes' entry containing 'param_summary' (per-component posterior mean/std/percentiles, same format as the bootstrap summary), 'samples' (equally weighted posterior draws for every parameter and derived quantity), 'logz', 'logzerr', 'noise' (posterior median noise standard deviation), and 'curves' (posterior percentile bands of the model).

---

## `unmix_coercivity_spectrum`

```python
unmix_coercivity_spectrum(x, spectrum, n_components=None, initial_parameters=None, vary_skew=DEFAULT_UNMIX_VARY_SKEW, weights=None, dp_bounds=(0.01, 2.0), skew_bounds=(-10.0, 10.0))
```

Unmix a coercivity spectrum into skew-normal (log-Gaussian) components.

Fits the derivative spectrum dM/dlog10(B) with a sum of skew-normal
densities, following the approach popularized by Kruiver et al. (2001)
and the MAX UnMix program (Maxbauer et al., 2016). Compared to fitting
the measured curve directly (unmix_backfield_curve) this operates on a
numerically differentiated (and possibly smoothed) version of the data,
so the choice of smoothing can influence the result; the advantage is
that components are fit in the space where they are most readily
interpreted visually.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT), e.g. midpoints from coercivity_spectrum_from_curve.
- **spectrum** (`array - like`) — Coercivity spectrum values at x (magnetization per decade).
- **n_components** (`int`) — Number of components. Required if initial_parameters is None.
- **initial_parameters** (`pandas.DataFrame`) — Initial guesses with columns 'contribution', 'location', 'dp', 'skew' (one row per component). If None, automatic estimates from estimate_coercivity_components are used.
- **vary_skew** (`bool`) — If False, skew values are fixed at their initial values (default 0, i.e. symmetric log-Gaussian components).
- **weights** (`array - like`) — Multiplicative weights applied to the residuals.
- **dp_bounds** (`tuple`) — (min, max) bounds on the dp scale parameter in log10 units.
- **skew_bounds** (`tuple`) — (min, max) bounds on the skew shape parameter.

**Returns**

- `dict` — Standardized result dictionary with keys including 'params' (a DataFrame of fitted parameters, linearized standard errors, and derived quantities such as B_mean_mT and proportion), 'y_fit', 'residuals', 'stats' (rss, r_squared, aic, bic, ...), 'success', and 'initial_parameters'.

---

## `unmixing_bootstrap`

```python
unmixing_bootstrap(result, n_boot=500, resample='cases', proportion=1.0, noise_level=None, random_seed=None, n_grid=200, verbose=False)
```

Bootstrap uncertainty estimation for an unmixing result.

Repeatedly refits the model to resampled data, starting each fit from
the best-fit parameters, and summarizes the distributions of parameters
and model curves. Two resampling schemes are available:

- 'cases': data points are drawn with replacement ('proportion'
  controls the resample size relative to the data). With
  proportion=0.95 and noise_level=0.02 this emulates the resampling
  scheme of the MAX UnMix program (Maxbauer et al., 2016).
- 'residuals': the best-fit curve is perturbed with resampled fit
  residuals, preserving the field spacing of the original data.

Bootstrap distributions capture the full nonlinearity of the model and
are generally more trustworthy than the linearized standard errors in
the 'params' table, especially for strongly overlapping components.

**Parameters**

- **result** (`dict`) — Result from unmix_coercivity_spectrum or unmix_backfield_curve.
- **n_boot** (`int`) — Number of bootstrap replicates (default 500).
- **resample** (`str`) — 'cases' or 'residuals'.
- **proportion** (`float`) — Fraction of the data resampled per replicate for 'cases' (default 1).
- **noise_level** (`float`) — If given, multiplicative Gaussian noise with this relative standard deviation is added to each resampled dataset (MAX UnMix uses 0.02).
- **random_seed** (`None, int, or numpy.random.Generator`) — Seed for reproducibility.
- **n_grid** (`int`) — Number of grid points for the model-curve confidence bands.
- **verbose** (`bool`) — Print a progress summary.

**Returns**

- `dict` — A copy of the input result with an added 'bootstrap' entry containing 'param_summary' (per-component mean/std/percentiles for each parameter and derived quantity), 'curves' (percentile bands of the total and per-component model on 'x_grid'), 'n_success', and 'param_samples' (the raw bootstrap parameter arrays).

---

## `unmixing_multistart`

```python
unmixing_multistart(x, magnetization, method=DEFAULT_UNMIX_METHOD, n_components=2, n_starts=100, vary_skew=DEFAULT_UNMIX_VARY_SKEW, curve_type='backfield', random_seed=None, location_tolerance=0.1, proportion_tolerance=0.05, verbose=False, kwargs={})
```

Map the distinct unmixing solutions reachable from many initializations.

Coercivity unmixing is a non-convex problem: fits started from
different initial parameters can converge to different local minima
that describe the data almost equally well. A single fit (and a
bootstrap of it, which restarts every replicate from the best-fit
values) is conditioned on one solution basin and therefore hides this
non-uniqueness. This function launches the fit from n_starts dispersed
random initializations (plus the automatic peak-detection estimate),
clusters the converged solutions, and reports each distinct solution
with its fit statistics and Akaike weight, making the degeneracy of
the decomposition explicit.

**Parameters**

- **x** (`array - like`) — log10 of field values (mT).
- **magnetization** (`array - like`) — Remanence curve values at x (e.g. 'magn_mass_shift').
- **method** (`str`) — Registered unmixing method used for each fit (default DEFAULT_UNMIX_METHOD, i.e. 'spectrum').
- **n_components** (`int`) — Number of components (default 2).
- **n_starts** (`int`) — Number of random initializations (default 100).
- **vary_skew** (`bool`) — Whether skew varies during fitting; random starts draw skew from [-5, 5] when True.
- **curve_type** (`str`) — 'backfield' or 'acquisition'.
- **random_seed** (`None, int, or numpy.random.Generator`) — Seed for reproducible starting points.
- **location_tolerance** (`float`) — Two solutions are considered the same when all component locations agree within this tolerance (log10 units, default 0.1) and all proportions agree within proportion_tolerance.
- **proportion_tolerance** (`float`) — Proportion agreement tolerance (default 0.05).
- **verbose** (`bool`) — Print a summary of the distinct solutions.
- ****kwargs** — Passed through to the unmixing method.

**Returns**

- `dict` — The best (lowest RSS) result dictionary, augmented with a 'multistart' entry containing 'solutions' (a DataFrame with one row per distinct solution: n_hits, rss, r_squared, aic, delta_aic, akaike_weight, and per-component B_mean_mT / sd_log / proportion columns), 'results' (the representative result dictionary for each solution, in the same order), 'n_starts', and 'n_converged'.

---

## `verwey_estimate`

```python
verwey_estimate(temps, mags, t_range_background_min=50, t_range_background_max=250, excluded_t_min=75, excluded_t_max=150, poly_deg=3, plot_zero_crossing=False, plot_title=None, measurement_marker='o', measurement_color='FireBrick', background_fit_marker='s', background_fit_color='Teal', magnetite_marker='d', magnetite_color='RoyalBlue', verwey_marker='*', verwey_color='Pink', verwey_size=10, markersize=3.5)
```

Estimate the Verwey transition temperature and remanence loss of magnetite from MPMS data.
Plots the magnetization data, background fit, and resulting magnetite curve, and 
optionally the zero-crossing.

**Parameters**

- **temps** (`pd.Series`) — Series representing the temperatures at which magnetization measurements were taken.
- **mags** (`pd.Series`) — Series representing the magnetization measurements.
- **t_range_background_min** (`int or float`) — Minimum temperature for the background fitting range. Default is 50.
- **t_range_background_max** (`int or float`) — Maximum temperature for the background fitting range. Default is 250.
- **excluded_t_min** (`int or float`) — Minimum temperature to exclude from the background fitting range. Default is 75.
- **excluded_t_max** (`int or float`) — Maximum temperature to exclude from the background fitting range. Default is 150.
- **poly_deg** (`int`) — Degree of the polynomial for background fitting. Default is 3.
- **plot_zero_crossing** (`bool`) — If True, plots the zero-crossing of the second derivative. Default is False.
- **plot_title** (`str`) — Title for the plot. Default is None.
- **measurement_marker** (`str`) — Marker symbol for measurement data. Default is 'o'.
- **measurement_color** (`str`) — Color for measurement data. Default is 'black'.
- **background_fit_marker** (`str`) — Marker symbol for background fit data. Default is 's'.
- **background_fit_color** (`str`) — Color for background fit data. Default is 'C1'.
- **magnetite_marker** (`str`) — Marker symbol for magnetite data. Default is 'd'.
- **magnetite_color** (`str`) — Color for magnetite data. Default is 'C0'.
- **verwey_marker** (`str`) — Marker symbol used to denote the Verwey transition estimate on the plot. Default is '*'.
- **verwey_color** (`str`) — Color of the marker representing the Verwey transition estimate. Default is 'Pink'.
- **verwey_size** (`int`) — Size of the marker used for the Verwey transition estimate. Default is 10.
- **markersize** (`float`) — Size of the markers. Default is 3.5.

**Returns**

- `float` — Estimated Verwey transition temperature.
- `float` — Estimated remanence loss. 

**Examples**

```python
>>> temps = pd.Series([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
>>> mags = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
>>> verwey_estimate(temps, mags)
(75.0, 0.5)
```

---

## `verwey_estimate_interactive`

```python
verwey_estimate_interactive(measurements, specimen, method, figsize=(11, 5))
```

Create an interactive widget for estimating the Verwey transition temperature from low temperature remanence measurements.

This function displays interactive sliders and controls for adjusting background fitting parameters
and temperature ranges, allowing the user to visually estimate the Verwey transition temperature (T_v)
for a selected specimen and measurement method. The function updates plots in real-time according to user input,
enabling exploration of parameter effects on the calculated transition.

**Parameters**

- **measurements** (`pandas.DataFrame`) — low temperature remanence measurement data containing temperature and magnetization columns for multiple specimens.
- **specimen** (`str or ipywidgets.Dropdown`) — Specimen to analyze, given either as a plain specimen name or as a selection widget (e.g. from ``verwey_specimen_method_selection_interactive``); for a widget the current ``.value`` is read when this function runs, so rerun the cell after changing the dropdown.
- **method** (`str or ipywidgets.Dropdown`) — Measurement method ('LP-FC' or 'LP-ZFC'), as a plain string or a selection widget.
- **figsize** (`tuple of (float, float)`) — Size of the matplotlib figure, by default (11, 5).

**Notes**

- The function uses `ipywidgets` for interactive controls and `matplotlib` for visualization.
- The background fit and excluded temperature ranges can be adjusted using sliders.
- The polynomial degree of the background fit is also adjustable.
- A reset button restores the default slider values.
- The function relies on supporting functions such as `extract_mpms_data_dc`, `thermomag_derivative`, and `calc_verwey_estimate`.

**Returns**

- `None` — This function is intended for use in Jupyter notebooks or environments that support interactive widgets and inline plotting. It displays interactive sliders and plots but does not return a value. 

**Examples**

```python
>>> verwey_estimate_interactive(measurements_df, specimen_dropdown, method_dropdown)
>>> verwey_estimate_interactive(measurements_df, 'NED2-8c', 'LP-FC')
Displays an interactive interface for estimating the Verwey transition temperature.
```

---

## `verwey_estimate_multiple_specimens`

```python
verwey_estimate_multiple_specimens(specimens_with_params, measurements)
```

Analyze Verwey transitions for a list of specimens with unique parameters.

This function uses either field-cooled (FC) or zero-field cooled (ZFC) data depending on the 
method_codes provided in each specimen's parameters. If "LP-FC" is found in the colon-delimited 
method_codes, FC data is used; if "LP-ZFC" is found, ZFC data is used.

**Parameters**

- **specimens_with_params** (`list of dict`) — List of specimen dictionaries. Each dictionary should contain:     - 'specimen_name' : str       The name of the specimen.     - 'params' : dict       Dictionary containing:         - 't_range_background_min' : int or float         - 't_range_background_max' : int or float         - 'excluded_t_min' : int or float         - 'excluded_t_max' : int or float         - 'poly_deg' : int         - 'method_codes' : str           Colon-delimited string that must include either "LP-FC" or "LP-ZFC".
- **measurements** (`object`) — Measurements dataframe in MagIC format.

**Returns**

- `pd.DataFrame` — DataFrame containing the Verwey transition estimates and the input parameters for each specimen. Columns include:     - 'specimen'     - 'critical_temp'     - 'critical_temp_type'     - 'remanence_loss' plus the additional parameters from the input. 

**Raises**

- `ValueError` — If neither "LP-FC" nor "LP-ZFC" is found in the method_codes for a specimen.
- `Exception` — Propagates exceptions raised during data extraction or analysis.

---

## `verwey_specimen_method_selection_interactive`

```python
verwey_specimen_method_selection_interactive(measurements)
```

Creates and displays dropdown widgets for selecting a specimen and the corresponding
available method codes (specifically 'LP-FC' and 'LP-ZFC') from a given DataFrame of measurements.
This function filters the measurements to include only those with desired method codes,
dynamically updates the method dropdown based on the selected specimen, and organizes
the dropdowns vertically in the UI.

Parameters:
    measurements (pd.DataFrame): The DataFrame containing measurement data with columns
                                 'specimen' and 'method_codes'. It is expected to have
                                 at least these two columns where 'specimen' identifies
                                 the specimen name and 'method_codes' contains the method
                                 codes associated with each measurement.

Returns:
    tuple: A tuple containing the specimen dropdown widget (`ipywidgets.Dropdown`)
           and the method dropdown widget (`ipywidgets.Dropdown`). The specimen dropdown
           allows for the selection of a specimen, and the method dropdown updates to
           display only the methods available for the selected specimen. The initial
           selection in the specimen dropdown is set to the first specimen option.

Note:
    The method dropdown is initially populated based on the methods available for the
    first selected specimen. The available methods are specifically filtered for 'LP-FC'
    and 'LP-ZFC' codes.

---

## `zero_crossing`

```python
zero_crossing(dM_dT_temps, dM_dT, make_plot=False, xlim=None, verwey_marker='*', verwey_color='Pink', verwey_size=10)
```

Calculate the temperature at which the second derivative of magnetization with respect to 
temperature crosses zero. This value provides an estimate of the peak of the derivative 
curve that is more precise than the maximum value.

The function computes the second derivative of magnetization (dM/dT) with respect to 
temperature, identifies the nearest points around the maximum value of the derivative, 
and then calculates the temperature at which this second derivative crosses zero using 
linear interpolation.

Parameters:
    dM_dT_temps (pd.Series): A pandas Series representing temperatures corresponding to
                             the first derivation of magnetization with respect to temperature.
    dM_dT (pd.Series): A pandas Series representing the first derivative of 
                       magnetization with respect to temperature.
    make_plot (bool, optional): If True, a plot will be generated. Defaults to False.
    xlim (tuple, optional): A tuple specifying the x-axis limits for the plot. Defaults to None.
    verwey_marker : str, optional
        Marker symbol used to denote the Verwey transition estimate on the plot. Default is '*'.
    verwey_color : str, optional
        Color of the marker representing the Verwey transition estimate. Default is 'Pink'.
    verwey_size : int, optional
        Size of the marker used for the Verwey transition estimate. Default is 10.

Returns:
    float: The estimated temperature at which the second derivative of magnetization 
           with respect to temperature crosses zero.

Note:
    The function assumes that the input series `dM_dT_temps` and `dM_dT` are related to 
    each other and are of equal length.
