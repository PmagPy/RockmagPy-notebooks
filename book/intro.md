---
title: RockmagPy
---

# RockmagPy: tools for rock magnetic data analysis

RockmagPy is a module of the [PmagPy](https://pmagpy.github.io/PmagPy-docs/intro.html) tools written in Python that is focused on the analysis of rock magnetic data. The functions within the RockmagPy module facilitate the visualization and interpretation of the results from experiments such as hysteresis loops, low-temperature remenance experiments, and thermomagnetic curves. 

RockmagPy is designed to work with the [MagIC database](https://earthref.org/MagIC) and enables manipulation of downloaded MagIC data sets as well as preparation of new contributions for uploading to the MagIC database. Utilizing RockmagPy functions within Jupyter notebooks enables fully documented and nicely illustrated data analysis. Example notebooks have been developed to aid in the use of RockmagPy. These notebooks are rendered within this JupyterBook. These notebooks can also be downloaded from the [RockmagPy-notebooks repository](https://github.com/pmagpy/rockmagpy-notebooks) to be run locally on your own computer or run in the cloud on the [EarthRef JupyterHub](https://jupyterhub.earthref.org/).

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

Launch RockmagPy notebooks on Earthref's JupyterHub. Run the `rockmag_set_up.ipynb` notebook to start.
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
The versions that are linked to here are static (that is that they don't run the code).
Some of the interactive features that have been developed are not visible. 
They do give an overall sense of what functions are available within RockmagPy and explain
some of the associated methods.

### Code/data nuts and bolts

- [**RockmagPy set-up**](./rockmag_set_up.ipynb) This notebook can be imported into JupyterHub (like the EarthRef JupyterHub) and be used to import the latest code and notebooks.
- [**Rockmag MagIC data unpack**](./rockmag_data_unpack.ipynb) This notebook takes data from rock magnetic experiments that are in MagIC format and splits out the different tables. It also summarizes all of the specimens and experiments within the measurements table.

### VSM experiments

These notebooks support the analysis of experiments conducted on a vibrating sample magnetometer (VSM).

- [*Hysteresis Loop Analysis*](../hysteresis.ipynb) *not yet developed*
- [*Back Field Curves*](../backfield.ipynb) *not yet developed*
- [*IRM Curves*](../irm_acquisition.ipynb) *not yet developed*
- [*Thermomagnetic Curves (Ms-T)*](../thermomagnetic.ipynb) *not yet developed*

### MPMS data

These notebooks support the analysis of low-temperature experiments conducted on a Magnetic Property Measurement System (MPMS).

- [**MPMS Data Plotting**](../MPMS_plot.ipynb) Within this notebook, MPMS data within a measurements file can be plotted. This plotting includes the calculation and visualization of derivative curves.
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