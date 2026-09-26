# Thermomagnetic data analysis

These notebooks support the analysis of thermomagnetic data.

- [**Thermomagnetic Curves**](./high_T_susceptibility.ipynb) This notebook plots and
  analyzes data from high-temperature susceptibility experiments.
- [**Curie temperature estimation**](./curie_temperature_estimation.ipynb) This notebook
  demonstrates the Curie temperature estimation methods implemented in RockmagPy, their
  physical basis, and the caveats that attend each of them — for both strong-field
  magnetization *Mₛ(T)* curves and low-field susceptibility *χ-T* curves.
- [**Evaluation of Curie temperature code in PmagPy**](./curie_code_evaluation.md) A
  written evaluation of the current and legacy Curie temperature estimation code across
  PmagPy, with recommendations on what to use.

## A note on methods and their applicability

Curie temperature estimates are method-dependent, and the appropriate methods differ
between data types. Estimates from different methods applied to the same curve differ
systematically — by several degrees to several tens of degrees
([Lattard et al., 2006](https://doi.org/10.1029/2006JB004591)) — so the estimation
method (and any smoothing or fitting windows) is part of the result and should be
reported with it.

For **strong-field magnetization *Mₛ(T)* curves**, the Landau-theory analysis of
[Fabian et al. (2013)](https://doi.org/10.1029/2012GC004440) places the Curie
temperature at the *inflection point* of the in-field curve (the minimum of
*dM/dT*), a position that is effectively independent of the applied field. The
classical maximum-curvature (second-derivative maximum) and "two-tangent"
([Grommé et al., 1969](https://doi.org/10.1029/JB074i022p05277)) estimates coincide
with each other and lie systematically *above* the inflection-point *T꜀*, typically by
10–15 °C. RockmagPy provides all of these — `inflection`, `max_curvature`,
`two_tangent`, and a Landau equation-of-state fit (`landau`) with formal uncertainty —
so that estimates can be compared and the bias structure made explicit.

For **low-field susceptibility *χ-T* curves**, the two-tangent construction is *not*
physically justified, and RockmagPy flags it when it is requested on susceptibility
data. [Fabian et al. (2013)](https://doi.org/10.1029/2012GC004440) explain this well:

> "The importance of the difference between determining *T<sub>c</sub>* from *Mₛ(T)* and *χ-T* is pointed out by [Petrovský and Kapicka (2006)](https://doi.org/10.1029/2006JB004507), where methods to determine *T<sub>c</sub>* from measurements of the initial susceptibility are analyzed. They conclude that the two-tangent method is not suitable for *χ-T* and can considerably overestimate *T<sub>c</sub>*. The physical origin of *χ-T* close to *T<sub>c</sub>* is more challenging than that of *Mₛ(T)*, because a number of low-field effects are important for *χ-T*, but become negligible in the higher fields used to infer *Mₛ(T)*. The variation of *m* depends not only on the variation of *Mₛ(H,T)* with field *H*, it also contains a contribution from a rotation of the ordered moment with respect to an easy magnetization axis, and contributions from thermally activated switching of small independent – but already magnetically ordered – regions (e.g., SP particles). In large bulk material, domain-wall movement contributes to *χ-T* even slightly below *T<sub>c</sub>*. In nanoparticles, the inhomogeneity of *Mₛ* due to the different exchange coupling of inner and surface atoms is of additional importance."

And from [Petrovský and Kapicka (2006)](https://doi.org/10.1029/2006JB004507):

> "...susceptibility for *T > T<sub>c</sub>* and *T < T<sub>c</sub>* increases to infinity, and we have to use analytical formulas developed for susceptibility behavior above the Curie point. Here, due to the geometry of the susceptibility curve, the two-tangent method will always yield temperature above the inflection point, which is higher than the temperature at which the substance starts to obey the paramagnetic Curie-Weiss law. The resulting error in *T<sub>c</sub>* (or *T<sub>N</sub><sup>1</sup>*) can be on the order of several degrees to several tens of degrees. Therefore, in the case of temperature dependence of magnetic susceptibility, application of the two-tangent method is not justified."

Furthermore:

> "In the case of synthetic magnetite and hematite, with sharp Hopkinson peak, the difference between transition temperatures determined using the two-tangent method and Curie-Weiss paramagnetic law is in the order of some few degrees. In the case of samples with wide susceptibility maximum and gradual decrease, reflecting e.g., wide distribution of grain sizes, or in the case of substituted hematite, application of the two-tangent method to susceptibility curves overestimates the transition temperature by several tens of degrees."

For *χ-T* data, the derivative-based estimates on the descending step and the
inverse-susceptibility (Curie-Weiss) extrapolation advocated by
[Petrovský and Kapicka (2006)](https://doi.org/10.1029/2006JB004507) are the
appropriate quantitative approaches — with the caveats that the Hopkinson peak marks
blocking rather than the Curie temperature, and that the Curie-Weiss extrapolation
yields the paramagnetic Curie temperature *θ* ≥ *T꜀*. The
[Curie temperature estimation notebook](./curie_temperature_estimation.ipynb)
demonstrates each of these methods and their caveats on measured data.
