---
title: RockmagPy
---

# RockmagPy: tools for rock magnetic data analysis

RockmagPy is a module of the [PmagPy](https://pmagpy.github.io/PmagPy-docs/intro.html) tools written in Python that is focused on the analysis of rock magnetic data. The functions within the RockmagPy module facilitate the visualization and interpretation of the results from experiments such as hysteresis loops, low-temperature remenance experiments, and thermomagnetic curves. 

RockmagPy is designed to work with the MagIC database (https://earthref.org/MagIC) and enables manipulation of downloaded MagIC data sets as well as preparation of new contributions for uploading to the MagIC database. Utilizing RockmagPy functions within Jupyter notebooks enables fully documented and nicely illustrated data analysis. Example notebooks have been developed to aid in the use of RockmagPy. These notebooks are rendered within this JupyterBook.

## RockmagPy notebooks

### VSM experiments

These notebooks support the analysis of experiments conducted on a vibrating sample magnetometer (VSM).

- [Hysteresis Loop Analysis](../hysteresis.ipynb) *not yet developed*
- [Back Field Curves](../backfield.ipynb) *not yet developed*
- [IRM Curves](../irm_acquisition.ipynb) *not yet developed*
- [Thermomagnetic Curves (Ms-T)](../thermomagnetic.ipynb) *not yet developed*

### MPMS data

These notebooks support the analysis of low-temperature experiments conducted on a Magnetic Property Measurement System (MPMS).

- [MPMS Data Plotting](../MPMS_plot.ipynb)
- [Verwey Transition Fitting](../MPMS_verwey_fit.ipynb)

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

Launch RockmagPy notebooks on the Earthref JupyterHub. Run the rockmag_set_up.ipynb notebook to get started.
:::

:::{grid-item-card}
:link: https://www2.earthref.org/MagIC/
:class-header: bg-light

MagIC database ⚡
^^^

Develop new contributions and analyze existing datasets in the MagIC database.
:::

::::

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