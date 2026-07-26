---
title: RockmagPy
---

# RockmagPy: tools for rock magnetic data analysis

RockmagPy is a module of the [PmagPy](https://pmagpy.github.io/PmagPy-docs/intro.html) tools written in Python that is focused on the analysis of rock magnetic data. The functions within RockmagPy facilitate the visualization and interpretation of the results from experiments such as hysteresis loops, low-temperature remenance cycling, and thermomagnetic curves. 

RockmagPy is designed to work with the [MagIC database](https://earthref.org/MagIC). It enables manipulation of MagIC data sets as well as the preparation of contributions for uploading to MagIC. Example notebooks have been developed to aid in the use of RockmagPy and with the goal of supporting well-documented data analysis. These notebooks can be downloaded from the [RockmagPy notebooks repository](https://github.com/pmagpy/rockmagpy-notebooks) to be run locally on one's computer.

## Jupyter notebooks using RockmagPy

Notebooks have been developed to use RockmagPy to analyze data from rock magnetic experiments.
The notebook versions that are shown within this book are static (i.e. they don't run the code).
Some of the interactive features that have been developed are not visible as they would be in a running notebook. 
The rendered notebooks do give an overall sense of what functions are available within RockmagPy and explain associated methods. These notebooks can be run by downloading them from the [RockmagPy repository](https://github.com/pmagpy/rockmagpy-notebooks) and running them in a Python environment on your computer. The [**Rockmag MagIC data unpack**](../rockmag_data_unpack.ipynb) notebook illustrates how data from rock magnetic experiments that are in MagIC format can be imported from the database.

The notebooks are organized into the following groups:

- [**MPMS data analysis**](../MPMS_notebooks/MPMS_intro.md) — analysis of low-temperature remanence and AC susceptibility data from Magnetic Property Measurement System (MPMS) experiments, including Verwey transition and goethite component fitting.
- [**VSM data analysis**](../hysteresis_backfield_notebooks/VSM_intro.md) — processing of hysteresis loops and backfield curves from vibrating sample magnetometer (VSM) experiments, including coercivity spectra unmixing.
- [**Thermomagnetic data analysis**](../thermomagnetic_notebooks/thermomag_intro.md) — analysis of high-temperature susceptibility data, including Curie temperature estimation.
- [**AMS data analysis**](../anisotropy_notebooks/anisotropy_intro.md) — visualization and statistical analysis of anisotropy of magnetic susceptibility (AMS) data.
- [**Rock Magnetic Bestiary**](../RMB_notebooks/RMB_intro.md) — interactive visualization of reference datasets from well-characterized materials developed at the Institute for Rock Magnetism.

Additional background material supports the interpretation of these data:

- [**Low-temperature behavior**](background/low_temperature_behavior.md) — background on low-temperature magnetic behavior including the Verwey and Morin transitions, superparamagnetism, and paramagnetism.
- [**Resources**](resources.md) — resources for learning more about rock magnetism.
- [**rockmag API reference**](rockmag_api.md) — reference sheet for the functions in the `pmagpy.rockmag` module that are used throughout these notebooks.

## Additional resources

::::{grid} 1 1 2 3

:::{grid-item-card}
:link: https://pmagpy.github.io/PmagPy-docs/installation/PmagPy_install

Installation 🪴
^^^

Instructions for how to install PmagPy on your computer.

:::

:::{grid-item-card}
:link: https://pmagpy.github.io/Essentials-JupyterBook/

Essentials of Paleomagnetism 📖
^^^

Explore the interactive Essentials of Paleomagnetism textbook that includes rock magnetism fundamentals.
:::

:::{grid-item-card}
:link: https://www2.earthref.org/MagIC/

MagIC database ⚡
^^^

Develop new contributions and analyze existing datasets in the MagIC database.
:::

::::

# Acknowledgements

::::{grid} 2 2 2 2

:::{grid-item}

```{image} images/logos/NSF_logo.png
:class: m-auto
:width: 150px
```

:::

:::{grid-item}
Current development of RockmagPy is supported by the National Science Foundation through its [support for the MagIC database](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148719).
:::

::::
