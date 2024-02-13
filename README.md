# RockmagPy notebooks
This repository is for the development and demonstration of the `rockmag.py` function module. The ultimate goal is to integrate additional functionality for the analysis of rock magnetic experiments into the `pmagpy` package. To facilitate development, the initial rockmag.py module is in this repository along with notebooks that use the rockmag.py functions to visualize and analyze data.

🚧 **WARNING: This repository is currently under initial development and is not yet ready for research or instructional use.** 🚧

## Notebooks

### General notebooks

#### [rockmag_set_up](./rockmag_set_up.ipynb)

This notebook can be imported into JupyterHub (like the EarthRef JupyterHub) and be used to import the lastest code and notebooks.

#### [rockmag_data_unpack](./rockmag_data_unpack.ipynb)

This notebook takes data from rock magnetic experiments that are in MagIC format and splits out the different tables. It also summarizes all of the specimens and experiments within the measurements table.

### MPMS notebooks

#### [MPMS_plot](./MPMS_plot.ipynb)

Within this notebook, MPMS data within a measurements file can be plotted. This plotting includes the calculation and visualization of derivative curves.

#### [MPMS_verwey_fit](./MPMS_verwey_fit.ipynb)

In this notebook, Verwey temperature estimates can be made from remanence upon warming experiments conducted on an MPMS.

## Contributing

While the repository is not yet ready for general use, contributions, suggestions, and feedback are welcome. _Thank you for your interest in contributing to the `rockmag.py` function library and the broader `pmagpy` project!_ 🌍🧲
