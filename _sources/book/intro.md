---
title: RockmagPy
---

# RockmagPy: tools for rock magnetic data analysis

RockmagPy is a module of the [PmagPy](https://pmagpy.github.io/PmagPy-docs/intro.html) tools written in Python that is focused on the analysis of rock magnetic data. The functions within RockmagPy facilitate the visualization and interpretation of the results from experiments such as hysteresis loops, low-temperature remenance cycling, and thermomagnetic curves. 

RockmagPy is designed to work with the [MagIC database](https://earthref.org/MagIC). It enables manipulation of MagIC data sets as well as the preparation of contributions for uploading to MagIC. Example notebooks have been developed to aid in the use of RockmagPy and with the goal of supporting well-documented data analysis. These notebooks can be downloaded from the [RockmagPy notebooks repository](https://github.com/pmagpy/rockmagpy-notebooks) to be run locally on one's computer or run in the cloud on the [EarthRef JupyterHub](https://jupyterhub.earthref.org/).

::::{grid} 1 1 2 3
:class-container: text-center
:gutter: 3

:::{grid-item-card}
:link: https://pmagpy.github.io/PmagPy-docs/installation/PmagPy_install
:class-header: bg-light

Installation 🪴
^^^

Instructions for how to install PmagPy on your computer.

:::

:::{grid-item-card}
:link: https://jupyterhub.earthref.org/
:class-header: bg-light

Live notebooks 🚀
^^^

Launch RockmagPy notebooks in the cloud on Earthref's JupyterHub. Run `rockmag_set_up.ipynb` to start.
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

Notebooks have been developed to use RockmagPy to analyze data from rock magnetic experiments.
The notebook versions that are shown within this book are static (i.e. they don't run the code).
Some of the interactive features that have been developed are not visible as they would be in a running notebook. 
The rendered notebooks do give an overall sense of what functions are available within RockmagPy and explain associated methods. These notebooks can be run by either downloading them from the [RockmagPy repository](https://github.com/pmagpy/rockmagpy-notebooks) and running them in a Python environment on your computer or by running them in the cloud such as on the [EarthRef JupyterHub](https://jupyterhub.earthref.org/). The [**RockmagPy set-up**](../rockmag_set_up.ipynb) notebook is the place to get started with running notebooks on Jupyterhub. The [**Rockmag MagIC data unpack**](../rockmag_data_unpack.ipynb) notebook illustrates how data from rock magnetic experiments that are in MagIC format can be imported from the database.

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
