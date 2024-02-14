---
title: RockmagPy
---

# RockmagPy: tools for rock magnetic data analysis

RockmagPy is a module of the [PmagPy](https://pmagpy.github.io/PmagPy-docs/intro.html) tools written in Python that is focused on the analysis of rock magnetic data. The functions within RockmagPy facilitate the visualization and interpretation of the results from experiments such as hysteresis loops, low-temperature remenance experiments, and thermomagnetic curves. 

RockmagPy is designed to work with the [MagIC database](https://earthref.org/MagIC). It enables manipulation of MagIC data sets as well as preparation of contributions for uploading to MagIC. Example notebooks have been developed to aid in the use of RockmagPy and to enable well-documented data analysis. These notebooks can be downloaded from the [RockmagPy notebooks repository](https://github.com/pmagpy/rockmagpy-notebooks) to be run locally on one's computer or run in the cloud on the [EarthRef JupyterHub](https://jupyterhub.earthref.org/).

::::{grid} 1 1 2 3
:class-container: text-center
:gutter: 3

:::{grid-item-card}
:link: https://pmagpy.github.io/PmagPy-docs/installation/PmagPy_install.html
:link-type: doc
:class-header: bg-light

Installation 🪴
^^^

Instructions for how to install PmagPy for your use case.

:::

:::{grid-item-card}
:link: https://jupyterhub.earthref.org/
:class-header: bg-light

Live notebooks 🚀
^^^

Launch RockmagPy notebooks on Earthref's JupyterHub. Run `rockmag_set_up.ipynb` to start.
:::

:::{grid-item-card}
:link: https://www2.earthref.org/MagIC/
:class-header: bg-light

MagIC database ⚡
^^^

Develop new contributions and analyze existing datasets in the MagIC database.
:::

::::

## Jupyter notebooks using RockmagPy

This page links to available notebooks that have been developed to use RockmagPy to analyze data.
The notebook versions that are shown within this book are static (they don't run the code).
Some of the interactive features that have been developed are not visible as they would be in a running notebook. 
They do give an overall sense of what functions are available within RockmagPy and explain
some associated methods.

### Code/data nuts and bolts

- [**RockmagPy set-up**](./rockmag_set_up.ipynb) This notebook can be imported into JupyterHub (like the EarthRef JupyterHub) and be used to import the latest code and notebooks.
- [**Rockmag MagIC data unpack**](./rockmag_data_unpack.ipynb) This notebook takes data from rock magnetic experiments that are in MagIC format and splits out the different MagIC tables. It also summarizes all the specimens and associated experiments within the measurements table.

### VSM data analysis

These notebooks support the analysis of data generated during experiments conducted on a vibrating sample magnetometer (VSM).

- [*Hysteresis Loop Analysis*](../VSM_hysteresis.ipynb) *not yet developed*
- [*Back Field Curves*](../VSM_backfield.ipynb) *not yet developed*

### Kappabridge data analysis

- [*Thermomagnetic Curves*](../kappa_thermomagnetic.ipynb) *not yet developed*
- [*Anisotropy of Magnetic Susceptibility*](../kappa_anisotropy.ipynb) *not yet developed*

### MPMS data analysis

These notebooks support the analysis of low-temperature experimental data conducted on a Magnetic Property Measurement System (MPMS).

- [**MPMS Data Plotting (DC)**](../MPMS_plot_dc.ipynb) Within this notebook, MPMS experiments where remanence was measured upon warming and cooling. This plotting includes the calculation and visualization of derivative curves.
- [**MPMS Data Plotting (AC)**](../MPMS_plot_ac.ipynb) *not yet developed*
- [**Verwey Transition Fitting**](../MPMS_verwey_fit.ipynb) In this notebook, Verwey temperature estimates can be made from remanence upon warming experiments conducted on an MPMS.

# Acknowledgements

::::{grid} 2 2 2 2

:::{grid-item}
:columns: 4

```{image} images/logos/NSF_logo.png
:class: m-auto
:width: 150px
```

:::

:::{grid-item}
:columns: 7
Current development of RockmagPy is supported by the National Science Foundation through its [support for the MagIC database](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148719).
:::

::::
