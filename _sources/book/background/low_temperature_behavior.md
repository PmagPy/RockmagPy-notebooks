# Low temperature behavior

Low-temperature experiments involving both magnetic remanence and susceptibility can give significant insight into the magnetic mineralogy of natural samples. Characteristic low-temperature phase transitions in magnetic minerals including magnetite ([Verwey transition](verwey.md); Verwey 1939; Jackson & Moskowitz 2020), hematite ([Morin transition](morin.md); Morin 1950; Özdemir et al. 2008), and pyrrhotite ([Besnus transition](besnus.md); Besnus & Meyer 1964; Rochette et al. 1990)) alter magnetic remanence enabling identification of these minerals. Other minerals transition from being paramagnetic at room-temperature to being ferromagnetic (in the broad sense of the term) at low-temperature transition (such as [siderite](../minerals/siderite.md)) that has a Néel temperature of 37 K). Low-temperature behavior also provides insight into fundamental characteristics of magnetic mineral assemblages, including grain size and domain state. Such changes include the fact that nanoparticles that are superparamagnetic at room temperature (300 K) can become stable single-domain grains at low temperatures, which can enable the quantification of grain size (Worm & Jackson 1999). In addition to stepwise changes in anisotropy associated with phase transitions, more progressive changes in anisotropy will manifest in magnetic properties, particularly for multidomain grains (Moskowitz, et al. 1998).

An instrument that is widely used in rock magnetism, condensed matter physics, and materials science for such experiments is a Magnetic Properties Measurement System (MPMS). This high-sensitivity Superconducting Quantum Interference Device (SQUID) magnetometer is designed to measure the magnetic properties of materials across a wide range of temperatures and magnetic fields.

## Low-temperature remanence experiments

One typical low-temperature remanence experiment used in rock magnetic investigations involves imparting a saturating isothermal remanent magnetization at room temperature (300 K; RTSIRM) to a specimen, then cycling it from room temperature to low temperature and back again (typically 300 K → 10 K → 300 K). In the MagIC data model, the cooling of an RTSIRM to low temperature is denoted as the method code `LP-CW-SIRM:LP-MC`, with subsequent warming denoted as `LP-CW-SIRM:LP-MW`. 

Another common type of experiment is one in which a sample is cooled in zero field to low temperature (e.g., 10 K) and then given a saturation isothermal remanent magnetization (ZFC-LTSIRM). Alternatively, a field can be maintained throughout the cooling (field-cooled; FC-LTSIRM). After either a ZFC-LTSIRM or an FC-LTSIRM is imparted, the specimen is then warmed in the absence of an applied field from the low-temperature state up to room temperature (typically 10 K → 300 K). In the MagIC data model, a remanence-upon-warming experiment following an FC-LTSIRM is denoted with the method code `LP-FC`, and that following a ZFC-LTSIRM with `LP-ZFC`.

The interpretation of low-temperature remanence experiments typically begins with visual inspection of plots of magnetization versus temperature and the associated first derivative. In *rockmagpy*, data associated with such experiments are extracted from MagIC-format measurement files and then visualized using static plots (via matplotlib) or interactive plots (via plotly). An example of this type of data visualization can be found within the [MPMS_plot_dc.ipynb](https://pmagpy.github.io/RockmagPy-notebooks/MPMS_notebooks/MPMS_plot_dc.html) notebook.

---

## References

Besnus, M. J., & Meyer, A. J. (1964). Magnétisme de la pyrrhotine naturelle. Nouvelles données expérimentales. In Proc. Int. Conf. Mag. (Vol. 20, pp. 507–511). [DOI not available]

Jackson, M. J., & Moskowitz, B. (2020). On the distribution of Verwey transition temperatures in natural magnetites. Geophysical Journal International, 224(2), 1314–1325. [doi:10.1093/gji/ggaa516](http://dx.doi.org/10.1093/gji/ggaa516)

Liu, Q., Banerjee, S. K., Jackson, M. J., Chen, F., Pan, Y., & Zhu, R. (2004). Determining the climatic boundary between the Chinese loess and palaeosol: Evidence from aeolian coarse-grained magnetite. Geophysical Journal International, 156(2), 267–274. [doi:10.1111/j.1365-246X.2003.02148.x](http://dx.doi.org/10.1111/j.1365-246X.2003.02148.x)

Morin, F. J. (1950). Magnetic susceptibility of α-Fe₂O₃ and α-Fe₂O₃ with added titanium. Physical Review, 78(6), 819–820. [doi:10.1103/PhysRev.78.819.2](http://dx.doi.org/10.1103/PhysRev.78.819.2)

Moskowitz, B., Jackson, M., & Kissel, C. (1998). Low-temperature magnetic behavior of titanomagnetites. Earth and Planetary Science Letters, 157(3–4), 141–149. [doi:10.1016/S0012-821X(98)00033-8](http://dx.doi.org/10.1016/S0012-821X(98)00033-8)

Özdemir, Ö., Dunlop, D. J., & Berquó, T. S. (2008). Morin transition in hematite: Size dependence and thermal hysteresis. Geochemistry, Geophysics, Geosystems, 9(10), Q10Z01. [doi:10.1029/2008GC002110](http://dx.doi.org/10.1029/2008GC002110)

Rochette, P., Fillion, G., Mattéi, J.-L., & Dekkers, M. J. (1990). Magnetic transition at 30–34 Kelvin in pyrrhotite: Insight into a widespread occurrence of this mineral in rocks. Earth and Planetary Science Letters, 98(3–4), 319–328. [doi:10.1016/0012-821X(90)90034-U](http://dx.doi.org/10.1016/0012-821X(90)90034-U)

Worm, H.-U., & Jackson, M. (1999). The superparamagnetism of Yucca Mountain Tuff. Journal of Geophysical Research: Solid Earth, 104(B11), 25415–25425. [doi:10.1029/1999JB900285](http://dx.doi.org/10.1029/1999JB900285)