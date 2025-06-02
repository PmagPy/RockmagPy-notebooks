- [**Thermomagnetic Curves**](../high_T_susceptibility.ipynb) This notebook plots and
  analyzes data from high-temperature susceptibility experiments.

**Nota bene:** We do not provide Curie Temperature estimation by the "two-tangent"
method, as proposed for spontaneous magnetization *Mₛ(T)* curves by
[Grommé et al. (1969)](https://doi.org/10.1029/JB074i022p05277).
[Fabian et al. (2013)](https://doi.org/10.1029/2012GC004440) explain this well:

> "The importance of the difference between determining *T<sub>c</sub>* from *Mₛ(T)* and
> *χ-T* is pointed out by [Petrovský and Kapicka (2006)](https://doi.org/10.1029/2006JB004507),
> where methods to determine *T<sub>c</sub>* from measurements of the initial
> susceptibility are analyzed. They conclude that the two-tangent method is not
> suitable for *χ-T* and can considerably overestimate *T<sub>c</sub>*. The physical
> origin of *χ-T* close to *T<sub>c</sub>* is more challenging than that of *Mₛ(T)*,
> because a number of low-field effects are important for *χ-T*, but become
> negligible in the higher fields used to infer *Mₛ(T)*. The variation of *m* depends
> not only on the variation of *Mₛ(H,T)* with field *H*, it also contains a
> contribution from a rotation of the ordered moment with respect to an easy
> magnetization axis, and contributions from thermally activated switching of small
> independent – but already magnetically ordered – regions (e.g., SP particles). In
> large bulk material, domain-wall movement contributes to *χ-T* even slightly below
> *T<sub>c</sub>*. In nanoparticles, the inhomogeneity of *Mₛ* due to the different
> exchange coupling of inner and surface atoms is of additional importance."

And from [Petrovský and Kapicka (2006)](https://doi.org/10.1029/2006JB004507):

> "...susceptibility for *T > T<sub>c</sub>* and *T < T<sub>c</sub>* increases to infinity,
> and we have to use analytical formulas developed for susceptibility behavior above
> the Curie point. Here, due to the geometry of the susceptibility curve, the
> two-tangent method will always yield temperature above the inflection point, which
> is higher than the temperature at which the substance starts to obey the paramagnetic
> Curie-Weiss law. The resulting error in *T<sub>c</sub>* (or *T<sub>N</sub><sup>1</sup>*)
> can be on the order of several degrees to several tens of degrees. Therefore, in the
> case of temperature dependence of magnetic susceptibility, application of the
> two-tangent method is not justified."

Furthermore:

> "In the case of synthetic magnetite and hematite, with sharp Hopkinson peak, the
> difference between transition temperatures determined using the two-tangent method
> and Curie-Weiss paramagnetic law is in the order of some few degrees. In the case of
> samples with wide susceptibility maximum and gradual decrease, reflecting e.g., wide
> distribution of grain sizes, or in the case of substituted hematite, application of
> the two-tangent method to susceptibility curves overestimates the transition
> temperature by several tens of degrees."