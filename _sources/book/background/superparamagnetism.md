## Superparamagnetism

```{admonition} Source of text
:class: note
This text on superparamagnetic behavior at low temperature is from the Fall 2009 issue of the IRM Quarterly: *Bowles, Julie; Jackson, Mike; Chen, Amy; Solheid, Peat. (2009). IRM Quarterly, Volume 19, Number 3 (Fall 2009). Cover article: Interpretation of Low-Temperature Data Part 1: Superparamagnetism and Paramagnetism. Retrieved from the University Digital Conservancy, https://hdl.handle.net/11299/171316.*
 with minor edits by Nick Swanson-Hysell.
```

Superparamagnetism (SP) describes the state of a single-domain-sized grain when thermal energy is sufficient to overcome barriers to a reversal of magnetization. Here the term “reversal” assumes uniaxial anisotropy with two minimum energy states having anti-parallel moment orientations. Energy barriers to magnetization changing arise from magnetocrystalline, magnetoelastic and/or shape anisotropy, all of which are proportional to grain volume (V). When the energy barriers are large with respect to thermal energy, the magnetization is “blocked” and the probability of spontaneous reversal approaches nil. But when the barriers are relatively low, thermal excitations can result in reversal of the magnetization over very short time scales, and the grain is in a superparamagnetic state. At a given temperature the volume at which a particle goes from being unblocked to blocked is known as the blocking volume (V<sub>b</sub>). For a given volume, we can block the grain by lowering the temperature (i.e. decreasing the available thermal energy) below the blocking temperature (T<sub>b</sub>).

```{figure} https://upload.wikimedia.org/wikipedia/commons/thumb/d/d8/Ferrofluid_close.jpg/800px-Ferrofluid_close.jpg
:name: ferrofluid_figure
:align: center
:width: 60%

Ferrofluid (containing SP iron particles) placed over a
rare earth magnet. Source: [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Ferrofluid_close.jpg).
```

More formally, we can think about SP blocking in terms of a particle’s characteristic relaxation time ($\tau$) and how rapidly a particle or assemblage of particles may approach equilibrium. Based on Néel theory of thermally activated magnetization (Néel, 1949) and the terminology of Dunlop and Özdemir (1997), the equilibrium magnetization ($M$) for aligned, uniaxial grains is given by:

```{math}
:label: eq_equilibrium_magnetization
M_{\mathrm{eq}}(H_{0}, T) = M_{s}\tanh(\alpha)
```

which is $\approx M_{s}\alpha$ for small $\alpha$ where,

$
\alpha = \frac{\mu_{0}\,V\,M_{s}\,H_{0}}{k_{B}\,T}
$

and the relaxation time ($\tau$) in weak or zero field is:
 
```{math}
:label: eq_relaxation_time
\tau = \tau_{0} \exp\!\Bigl(\frac{K V}{k_{B} T}\Bigr)
```

where $M_{s}$ is the saturation magnetization, $\mu_{0}$ is the permeability of free space, $H_{0}$ is a (weak) applied field, and $\tau_{0}$ is a time constant related to the atomic reorganization time (on the order of $10^{-9}$ s), $k_{B}$ is Boltzmann’s constant, and $K$ is an anisotropy constant also related to particle microcoercivity. In the equilibrium state, some fraction of the population has the high‐energy moment orientation (antiparallel to the applied field) while others occupy the low‐energy orientation, and spontaneous reversals occur at the same rate in each orientation (like a two‐way reaction in chemical equilibrium). For an assemblage of randomly oriented single‐domain grains, the relationship follows a Langevin function as:

```{math}
:label: eq_sd_Langevin
M_{\mathrm{eq}}(H_{0}, T) \;=\; M_{s}\,\Bigl[\coth(\alpha) \;-\; \frac{1}{\alpha}\Bigr]
```

$\approx \frac{M_s \alpha}{3}$ for small $\alpha$.

When $\tau$ is short compared to the observation time (the superparamagnetic state), the sample moment can equilibrate quickly with changes in applied field. In contrast, when $\tau$ is long the sample is blocked and in a stable single domain (SSD) state; it can take a very long time (millions or billions of years) to reach magnetic equilibrium. Whether or not a particle is SP or SSD will depend on both temperature and volume, through their effect on $\tau$. 

```{seealso}
A much more thorough theoretical discussion of blocking behavior can be found in [Dunlop and Özdemir (1997)](../resources.md).
```

## Remanence Data

The typical remanence experiments performed by researchers at the IRM include low-temperature cycling (300K → 20K → 300K) of a room-temperature (300K) saturation isothermal remanence (SIRM). A second common experiment cools the sample in zero field (ZFC) to low (e.g. 10K) temperature, where an SIRM is imparted. This low-temperature remanence is then measured as the sample warms back to room temperature. A room-temperature SIRM is not typically useful for studying SP particles, because they will not block at room temperature. However, many grains that are SP at room temperature will be SSD (blocked) at 20K. A sample with a low-temperature SIRM that warms in zero field ($M_{eq} = 0$) will demagnetize via the unblocking of these particles as they pass through their blocking temperature. If the sample contains a relatively narrow grain-size distribution of SP particles, the low-temperature remanence will unblock over a relatively narrow temperature interval (e.g. {numref}`SP_ZFC`, circles). If, however, the sample has a more distributed grain-size distribution, the unblocking will be more gradual (e.g. {numref}`SP_ZFC`, squares and diamonds). See Worm and Jackson (1999) for a discussion of volume distribution calculations based on the unblocking of a low-temperature remanence.

```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/SP_ZFC.png
:name: SP_ZFC
:align: center
:width: 60%

Saturation IRM acquired at 10K (ZFC) and measured on warming from 10K to 300K. The decrease in magnetization results from thermal unblocking of an SP population. Blue circles are from a Tiva Canyon Tuff sample with a relatively narrow grain-size distribution. Data from the same sample are shown in {numref}`SP_induced` and {numref}`SP_ac`b. Red squares and green triangles are from synthetic, glassy basalts with more distributed grain-size distributions.
```

However, not every such gradual decrease of low-temperature temperature remanence on warming is related to superparamagnetism. Moskowitz et al. (1998) demonstrate that multi-domain, high-Ti titanomagnetite exhibits similar behavior resulting from domain reorganization as magnetocrystalline anisotropy decreases on warming.

## Induced Magnetization Data

```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/SP_induced.png
:name: SP_induced
:align: center
:width: 50%

Induced FC-ZFC experiment on Tiva Canyon Tuff sample.
Sample is cooled in zero field and then measured on warming in a 5 mT DC field (blue circles). Sample is then cooled in a 5 mT field and measured on warming in a 5 mT field. Sample is fully unblocked when the two curves coincide. Data from the same sample are shown in {numref}`SP_ZFC` (circles) and {numref}`SP_ac`b.
```

In an experiment more common to the physics community (but also used by rock magnetists, e.g. Berquó et al., 2007), a sample is cooled in zero field (ZFC) to a low temperature such as 20K, and the induced sample moment ($M_i$) is then measured on warming in the presence of a small (e.g. 5 mT) field. This is followed by cooling the sample in the same 5 mT field (FC) and again measuring $M_i$ on warming in a 5 mT field ({numref}`SP_induced`). In the FC case, for \(T < T_b\), magnetization can be described by regular TRM theory, i.e. 

```{math}
:label: eq_TRM_theory
M_{TRM}(T) = M_{eq}(T_B) \frac{M_S(T)}{M_S(T_B)}
```
 
for a uniform distribution of particles with a single $T_B$. For non-uniform distributions (e.g. {numref}`SP_induced`), Eq. {eq}`eq_TRM_theory` must be integrated over all $T_B$. In the ZFC case, as the sample warms back up, $M_{i(ZFC)}$ is initially less than $M_{i(FC)}$ because the sample has blocked in a zero-field remanence. As the sample warms through $T_B$, $\tau$ becomes short enough that the moment approaches $M_{eq}$ in the applied field. The point where the FC and ZFC curves join together defines the (maximum) blocking temperature, and for $T > T_B$, both curves represent $M_{eq}$. If the sample has a narrow-enough grain-size distribution and enough is known about its properties in order to estimate $K$, particle volume can be estimated by solving Eq. {eq}`eq_sd_Langevin` for $V$.  

## AC Susceptibility Data  

In measuring low-field susceptibility ($\chi$), a weak, time-varying alternating field [$H(t) = H_0 \cos(\omega t)$] is applied to the sample, and the magnetic response is measured. Consider the behavior of a superparamagnetic grain above and below its blocking temperature.  At $T_B << T < T_C$, the grain is ferromagnetically ordered, but $\tau$ is very short compared to the observation time ($t_{obs}$), and the moment will respond exactly in phase with the alternating field.  

For an assemblage of randomly-oriented identical grains, susceptibility falls off with 1/T according to:  

```{math}
:label: eq_sp_chi
\chi_{sp} = \frac{\mu_0 V M_S^2}{3kT} = \frac{M_{eq}(H_0)}{H_0}
```

This type of behavior can be observed in {numref}`SP_ac`a at $T > \sim200K$.  

At $T << T_B$, $\tau$ is very long compared to $t_{obs}$, the grains are fully blocked in the SSD state, and susceptibility will be negligible in weak fields:  

```{math}
:label: eq_sd_chi
\chi_{sd} = \frac{2M_S}{3H_k}
```

where $H_k$ is the microcoercivity of the grain (the field required to reverse the magnetization in the absence of any thermal effects). This type of behavior is observed in {numref}`SP_ac`c at $T < \sim100K$.

```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/SP_ac.png
:name: SP_ac
:align: center
:width: 100%

AC susceptibility measured as a function of frequency and temperature on three samples of Tiva Canyon Tuff. Average grain size (determined by TEM) given above plots. As the dominant grain size increases, so does the temperature of the peak in both in-phase and out- of-phase susceptibility, which corresponds to the blocking temperature.
The sample shown in panel b with average grain size of 37 x 7.5 nm  is the same as in {numref}`SP_ZFC` (circles) and {numref}`SP_induced`.
```

In between the fully unblocked state (Eq. {eq}`eq_sp_chi`) and the fully blocked state (Eq. {eq}`eq_sd_chi`), susceptibility will be at a maximum near $T_B$ ({numref}`SP_ac`), where $\tau \approx t_{obs}$. Susceptibility will now depend both on $\tau$ and on the frequency ($f = \omega / 2\pi$) of the AC field. Also in this interval ($\tau \approx t_{obs}$), the sample moment will tend to lag behind the alternating field, and the measured signal can be broken down into a component that is in-phase with the applied field ($\chi'$) and a component that is out-of-phase ($\chi''$):

```{math}
:label: eq_chi_prime
\chi' = \frac{\chi_{sp}}{1 + \omega^2 \tau^2}
```

```{math}
:label: eq_chi_double_prime
\chi'' = \frac{\omega \tau \chi_{sp}}{1 + \omega^2 \tau^2}
```

and $(\chi’)^2 + (\chi’’)^2 = (\chi_{sp})^2 / (1 + \omega^2 \tau^2)$. At low $f$ and short $\tau$ ($\omega^2 \tau^2 \ll 1$) the moment can “keep up” with the driving field, and the out-of-phase signal will be small. At higher $f$ and longer $\tau$, the magnetization will lag farther behind, leading to a larger out-of-phase signal. In this way, we also see that both the in-phase and out-of-phase susceptibility components are frequency-dependent as the sample passes through the blocking temperature interval ({numref}`SP_ac`).

```{figure} https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/SP_ac.png
:name: SP_ac_basalt
:align: center
:width: 100%

AC susceptibility measured as a function of frequency
and temperature on same synthetic basaltic sample shown in {numref}`SP_ZFC`
(squares). Data show evidence for a broad SP grain-size distribution.
Additionally, the rapid decrease in susceptibility from 10K to 50K
results from a significant paramagnetic fraction.
```

As in the unblocking of a low-temperature SIRM, a discrete, narrow grain-size population will result in a relatively narrow temperature range over which frequency-dependence is observed and in a well-defined maximum in susceptibil-
ity (e.g. {numref}`SP_ac`a). This peak shifts to higher temperatures for a larger dominant grain size ({numref}`SP_ac`b,c). A broad grain-size distribution will result in a broad peak and/or frequency-dependence over a large temperature interval (e.g. {numref}`SP_ac_basalt`). 

```{seealso}
See Worm (1998) for forward calculations of expected AC susceptibility responses for a range of SP magnetite grain size distributions. See Shcherbakov and Fabian (2005) for a discussion of the inversion of AC susceptibility data for volume distributions, which the authors find to give a non-unique solution. Note also that we have thus far treated particle magnetic moment and microcoercivity as known quantities. In reality, however, they are known unknowns, especially for small particles with complicated surface spin arrangements arising from partial-oxidation, surface interactions, or reduced coordination of surface spins. See Egli (2009) for susceptibility-based grain volume distribution reconstruction that bypasses assumptions on particle magnetic moments and microcoercivity; also see “Thermal Fluctuation Tomography” below.
```

```{caution}
Frequency-dependence ($\chi_{fd}$) and the presence of out-of-phase susceptibility ($\chi''$) are often taken as indicators of a significant SP population. However, note that $\chi_{fd}$ can also arise from electrical eddy currents induced in the sample, and that $\chi''$ can be attributed to eddy currents or low-field hysteresis. See Jackson (2004; *IRM Quarterly* vol 13 no. 4) for a much more thorough discussion of out-of-phase (or “imaginary”) susceptibility. Also note that the absence of frequency-dependence at room temperature does not preclude the presence of an SP population; as shown in Fig. {numref}`SP_ac`, if $T_B$ is considerably below room temperature, $\chi_{fd}$ at room temperature will still be zero.
```

## Thermal fluctuation tomography

Another approach for jointly characterizing SP volume and microcoercivity ($H_K$) distributions is thermal fluctuation tomography (Jackson et al., 2006).

```{note}
This section is under construction and the reader is currently referred to the IRM Quarterly version of this article for it: https://hdl.handle.net/11299/171316 
```

## Superparamagnetic Imposters

The caution above noted other phenomena that might be mistaken for superparamagnetic behavior in the absence of additional information, and we briefly re-iterate this caution here. It is important to keep in mind that non-zero out-of-phase susceptibility can result from eddy currents or low-field hysteresis, in addition to thermally activated processes such as SP blocking or unblocking. One way to identify a thermally-activated process is to measure AC susceptibility at at least three frequencies; Néel theory predicts the relationship between $\chi_{fd}$ and $\chi''$ to be:

```{math}
:label: eq_chi_double_prime_neel
\chi'' = -\frac{\pi}{2} \frac{d\chi'}{d \ln(f)}
```

for non-interacting populations with smooth distributions of $\tau$ (see Shcherbakov and Fabian, 2005; Jackson, 2004). Note, however, that *any* thermally-activated process will obey this relationship; some electron-ordering phenomena in titanomagnetites, for example, appear to be thermally-activated (see, e.g. Carter-Stiglitz et al., 2006; Lagroix et al., 2004).

The gradual decrease in a low-temperature remanence on warming can be indicative of SP unblocking, but may also result from domain-wall reorganization (Moskowitz et al., 1998) or from partial oxidation. For the latter, the remanent magnetization decay <50 K, while still poorly understood, is probably related to additional anisotropy from surface-core coupling. The phase transition from spin glass to a ferri- or para-magnetic state often observed in titano-hematites is also accompanied by a dramatic decrease in low-temperature remanence on warming at temperatures <50–60 K. This transition is additionally characterized by frequency-dependent peaks in both $\chi'$ and $\chi''$ (Burton et al., 2008), but which do not obey Eq. {eq}`eq_chi_double_prime_neel`.

Finally, samples with a large paramagnetic/ferromagnetic ratio may exhibit a similar decrease resulting not from any real change in remanence, but from variation in an induced paramagnetic signal in the presence of an imperfectly-zeroed field (see the [paramagnetic low-temperature behavior section](./paramagnetism.md)).

## References

Berquó, T.S., Banerjee, S.K., Ford, R.G., Penn, R.L., & Pichler, T. (2007). High-crystallinity Si-ferrihydrite: An insight into its Néel temperature and size dependence of magnetic properties. *J. Geophys. Res.*, 112, B02S02. [doi:10.1029/2006JB004264](http://dx.doi.org/10.1029/2006JB004264)

Burton, B.P., Robinson, P., McEnroe, S.A., Fabian, K., & Boffa Ballaran, T. (2008). A low-temperature phase diagram for ilmenite-rich compositions in the system Fe₂O₃–FeTiO₃. *American Mineralogist*, 93, 870–872. [doi:10.2138/am.2008.2804](http://dx.doi.org/10.2138/am.2008.2804)

Carter-Stiglitz, B., Moskowitz, B., Solheid, P., Berquó, T.S., Jackson, M., & Kosterov, A. (2006). Low-temperature magnetic behavior of multidomain titanomagnetites: TM0, TM16, and TM35. *J. Geophys. Res.*, 111, B12S05. [doi:10.1029/2006JB004561](http://dx.doi.org/10.1029/2006JB004561)

Dunlop, D.J. (1965). Grain distributions in rocks containing single-domain grains. *J. Geomag. Geoelec.*, 17, 459–471. [doi:10.5636/jgg.17.459](http://dx.doi.org/10.5636/jgg.17.459)

Dunlop, D.J., & Özdemir, Ö. (1997). *Rock Magnetism: Fundamentals and Frontiers*. Cambridge University Press, New York.

Egli, R. (2009). Magnetic susceptibility measurements as a function of temperature and frequency I: inversion theory. *Geophys. J. Int.*, 177, 395–420. [doi:10.1111/j.1365-246X.2009.04081.x](https://doi.org/10.1111/j.1365-246X.2009.04081.x)

Jackson, M. (2004). Imaginary Susceptibility: A Primer. *IRM Quarterly*, 13(4). [https://conservancy.umn.edu/handle/11299/169889](https://conservancy.umn.edu/handle/11299/169889)

Jackson, M., Carter-Stiglitz, B., Egli, R., & Solheid, P. (2008). Characterizing the superparamagnetic grain distribution f(V,Hk) by thermal fluctuation tomography. *J. Geophys. Res.*, 113, B12S07. [doi:10.1029/2008JB005958](https://doi.org/10.1029/2008JB005958)

Lagroix, F., Banerjee, S.K., & Jackson, M.J. (2004). Magnetic properties of the Old Crow tephra: Identification of a complex iron titanium oxide mineralogy. *J. Geophys. Res.*, 109, B01104. [doi:10.1029/2003JB002678](https://doi.org/10.1029/2003JB002678)

Moskowitz, B.M., Jackson, M., & Kissel, C. (1998). Low-temperature magnetic behavior of titanomagnetites. *Earth Planet. Sci. Lett.*, 157, 141–149. [doi:10.1016/S0012-821X(98)00033-8](https://doi.org/10.1016/S0012-821X(98)00033-8)

Néel, L. (1949). Théorie du trainage magnétique des ferromagnétiques en grains fins avec application aux terres cuites. *Ann. Géophys.*, 5, 99–136.

Schlinger, C.M., Veblen, P.R., & Rosenbaum, J.G. (1991). Magnetism and magnetic mineralogy of ash flow tuffs from Yucca Mountain, Nevada. *J. Geophys. Res.*, 96, 6035–6052. [doi:10.1029/90JB02653](https://doi.org/10.1029/90JB02653)

Shcherbakov, V., & Fabian, K. (2005). On the determination of magnetic grain-size distributions of superparamagnetic particle ensembles using the frequency dependence of susceptibility at different temperatures. *Geophys. J. Int.*, 162, 736–748. [doi:10.1111/j.1365-246X.2005.02603.x](https://doi.org/10.1111/j.1365-246X.2005.02603.x)

Worm, H.-U., & Jackson, M.J. (1999). The superparamagnetism of the Yucca Mountain Tuff. *J. Geophys. Res.*, 104, 25,415–25,425. [doi:10.1029/1999JB900285](https://doi.org/10.1029/1999JB900285)

Worm, H.-U. (1998). On the superparamagnetic-stable single domain transition for magnetite, and frequency dependence of susceptibility. *Geophys. J. Int.*, 133, 201–210. [doi:10.1046/j.1365-246X.1998.1331468.x](https://doi.org/10.1046/j.1365-246X.1998.1331468.x)