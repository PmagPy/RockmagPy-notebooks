# Low-Temperature Magnetic Behavior of Siderite (FeCO$_3$)

## Background and Importance  
Siderite (iron(II) carbonate, FeCO$_3$) is a common diagenetic mineral in reducing environments that forms an important sink for iron {cite}`Dekkers2023a`. There can be substitution of Mn and Mg for Fe such that there is a solid solution series to rhodochrosite (MnCO$_3$) and magnesite (MgCO$_3$). It can be difficult to detect siderite with methods such as X-ray diffraction (XRD) in natural samples given fine grain sizes and low concentrations {cite}`Dekkers2023`. A sensitive way to detect siderite in sediments and sedimentary rocks is through its low-temperature magnetic properties {cite}`Housen1996a`. Siderite is paramagnetic at room-temperature, but it undergoes a **Néel transition** (antiferromagnetic ordering) at a low critical temperature $T_\mathrm{N}$ of ~37 K {cite}`Frederichs2003a`. Rhodochrosite has a similar paramagnetic to antiferromagnetic transition at ~34 K {cite}`Frederichs2003a`.

## The ~37 K Néel Temperature of Siderite  
**Pure siderite orders magnetically at $T_\mathrm{N}\approx37$ K**. This transition is one of a trio of low-temperature paramagnetic to antiferromagnetic transitions that occur in authigenic Fe-bearing minerals with the others being: rhodochrosite (MnCO$_3$) that orders at ~32–34 K and vivianite (Fe$_3$(PO$_4$)$_2\cdot8$H$_2$O) at ~12 K {cite}`Frederichs2003a`. In rocks, Fe$^{2+}$–Mn$^{2+}$ solid solutions can modify siderite's transition temperature between ~32 and 38 K.

```{figure} ../images/siderite_Roman-01_MPMS.png
:name: siderite_FC_ZFC_RTSIRM
:align: center
:width: 100%
Low-temperature remanence data developed on the specimen *siderite_Roman-01* from the Rock Magnetic Bestiary (as processed in this [notebook](https://pmagpy.github.io/RockmagPy-notebooks/RMB_notebooks/RMB_siderite.html)). There is a sharp remanence decay when warming during the field-cooled experiment (FC). This remanence loss is much larger than in the zero-field-cooled (ZFC) experiment. This difference is due to the additional thermal remanence acquired when the field is on as the specimen cools through the Néel temperature. In the RTSIRM experiment, there is change in magnetization while cooling below the Néel temperature. This behavoir is resulting from a small bias field within the MPMS instrument that is imparting magnetization. Note that this behavoir has the potential to be mistaken for the Besnus transition in pyrrhotite. However, the large difference between the FC and ZFC (FC/ZFC=17.4 for this specimen) is diagnostic of siderite and also helps to distinguish it from rhodochrosite.
```

```{iframe} ../images/siderite_Roman-01.html
:label: siderite_FC_ZFC_RTSIRM
:align: center
:width: 100%
:height: 400px
```

In the presence of an applied field, siderite and rhodochrosite acquire remanence below $T_\mathrm{N}$ due to spin-canted antiferromagnetism {cite}`Frederichs2003a`. If a sample is cooled in zero field to low temperature, given a saturation isothermal remanent magnetization (SIRM) at ~5–10 K, and then warmed (a zero-field-cooled or ZFC remanence experiment), the ~37 K transition causes a loss of magnetization ({numref}`siderite_FC_ZFC_RTSIRM`). The similarity between $T_\mathrm{N}$ of siderite and rhodochrosite can make it difficult to distinguish the two based on the transition temperature alone {cite}`Dekkers2023a`. If the sample is cooled from room-temperatures in the presence of a strong field (field-cooled, FC) before the remanence is measured the magnetization acquired (and then subsequently lost upon warming) is much larger than in a ZFC experiment ({numref}`siderite_FC_ZFC_RTSIRM`). This pronounced divergence between the FC and the ZFC curves (with FC>ZFC) is a hallmark of siderite ({numref}`siderite_FC_ZFC_RTSIRM`). This large difference is diagnostic for siderite relative to rhodochrosite and vivianite with a FC/ZFC ratio of 3 or greater proposed to identify siderite {cite}`Dekkers2023a`.

This FC to ZFC difference also plays an important role in distinguishing siderite’s transition from the Besnus transition of pyrrhotite that occurs at a similar temperature {cite}`Housen1996a`. Monoclinic pyrrhotite does not show a remanence divergence across in FC/ZFC data at its 32–34 K transition. Pyrrhotite’s low-temperature transition is a crystallographic/anisotropy change, but it remains a ferromagnetic (ferrimagnetic) mineral throughout, and experiments have found that ZFC vs. FC warming curves for pyrrhotite are nearly superposed in magnitude {cite}`Horng2018a`.

```{bibliography}
```