# Paramagnetism at low temperature

:::{note} Text source
This text on the paramagnetic at low temperature is from the Fall 2009 issue of the IRM Quarterly
Bowles, Julie; Jackson, Mike; Chen, Amy; Solheid, Peat. (2009). IRM Quarterly, Volume 19, Number 3 (Fall 2009). Cover article: Interpretation of Low-Temperature Data Part 1: Superparamagnetism and Paramagnetism. Retrieved from the University Digital Conservancy, https://hdl.handle.net/11299/171316.
:::

Any sample that has a paramagnetic phase in addition to ferromagnetic phases may result in a low-temperature signal that is the sum of both a ferromagnetic and a paramagnetic response. Because paramagnetic susceptibility (χₚ = M/H = C/T, where C is the Curie constant) is inversely proportional to temperature, it may dominate at very low temperatures when magnetization is blocked in other grains including superparamagnetic ones that therefore have very low susceptibility. Figure 1 shows AC‐susceptibility data for a sample with a high paramagnetic/ferromagnetic ratio, and at T < ~100 K, χₚ dominates with its characteristic 1/T signal.


```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/paramag_suscep_low_T.png
:label: paramag
:align: center
:width: 40%

AC susceptibility measured as a function of frequency and temperature on a synthetic basaltic sample. Data show evidence for a broad superparamagnetic grain-size distribution. Additionally, the rapid decrease in susceptibility from 10 K to 50 K results from a significant paramagnetic fraction.
```

In theory, there can be no paramagnetic contribution to remanence data because H = 0. However, in practice it is difficult to produce a true zero‐field environment inside the MPMS. It is not uncommon to have a “zero‐field” of ±1–2 µT, and for samples with a high paramagnetic/ferromagnetic ratio this is sufficient to induce a noticeable paramagnetic moment at low temperatures (Figure 2). This difficulty is exacerbated within the MPMS-3 instrument where a low-field is more difficult to obtain than the MPMS-XL.


```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/paramag_low_T.png
:label: paramag
:align: center
:width: 40%

Example of the effects of non-zero field in MPMS during remanence measurements. Panels show data from two different synthetic basaltic glass samples. A room-temperature SIRM was measured on cooling (red squares), and a 10 K SIRM was measured on warming (blue circles). At temperatures below ~50-100 K, the induced paramagnetic signal dominates, resulting from a small, residual field in the MPMS can be either negative (top), or positive (bottom).
```

