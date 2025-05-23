# VSM data analysis

These notebooks support the analysis of data generated during experiments conducted on a vibrating sample magnetometer (VSM). The coercivity spectra unmixing notebooks can also be applied to data from IRM acquisition or AF demagnetization experiments conducted on other magnetometers.

- [**Hysteresis Loop Analysis**](./hysteresis_processing.ipynb) This notebook processes raw hysteresis loop data following the protocols of {cite}`Jackson2010a`. Summary hysteresis parameters are calculated from a collection of loops.
- [**Hysteresis Loop Analysis (detailed work through)**](./hysteresis_processing_walk_through.ipynb) This notebook provides a step-by-step walk through of the processing of hysteresis loops. Each processing step is illustrated and described. This walk through can enable more granular data processing and illustrates the decision tree that is automatically implemented in [hysteresis loop analysis notebook](./hysteresis_processing.ipynb).
- [**Back Field Processing**](./backfield_processing.ipynb) This notebook analyses backfield data in order to obtain coercivity of remanence (B$_{cr}$) values. 
- [**Back Field Coercivity Spectra Unmixing**](./coercivity_unmixing_MaxUnmix.ipynb) This notebook implements the MaxUnmix algorithm for coercivity spectra analysis. 

```{bibliography}
```