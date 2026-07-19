# ---
# Appendix cells extracted from high_T_susceptibility.ipynb prior to merging the
# jupyterbook2 migration, which rewrote that notebook. Kept here so the local
# work isn't lost in the merge; re-insert into the updated notebook afterward.
#
# Format: jupytext "percent" — `# %% [markdown]` and `# %%` delimit cells, so
# this file can be reopened as a notebook (jupytext) or the cells copied back in.
#
# Depends on `selected_experiment` (and `rmag`) defined earlier in the notebook.
# ---

# %% [markdown]
# ## Appendix: determining the best temperature window for smoothing
#
# Smoothing of <span style="font-family: 'Times New Roman';">*$\chi$-T*</span> data prior to analysis is optional and is applied through the `smooth_window` parameter of `rmag.plot_X_T()` (a running average with the given temperature window; the default of `smooth_window=0` applies no smoothing). Rather than picking a window size arbitrarily, the window size can be chosen based on the trade-off between smoothness and misfit.
#
# The `rmag.optimize_moving_average_window()` function plots the **model roughness** against the **root mean square error** of the model **(note that the y-axis is inverted)** (a measure of the difference between the running average and the data). Each point on the graph is for a different temperature window size. The `knee` in the curve (the break of slope) can be interpreted as the optimal temperature window for smoothing. At this point (window) there is a decent trade off between the smoothness of the model (the running average) and error on the model.
#
# Set values for `min_temp_window` and `max_temp_window` to set the temperature range for exploring the rms and variance tradeoff.

# %%
min_temp_window = 0
max_temp_window = 50
steps = 20
fig, ax = rmag.optimize_moving_average_window(selected_experiment, min_temp_window, max_temp_window, steps=steps)

# %% [markdown]
# Upon inspecting the plot, a smoothing window size of about 10 K seems to be a good choice for this dataset. The chosen window can then be applied by passing it to the `smooth_window` parameter of `rmag.plot_X_T()`:

# %%
rmag.plot_X_T(selected_experiment,
              temp_unit='C',
              remove_holder=True,
              plot_derivative=True,
              smooth_window=10)
