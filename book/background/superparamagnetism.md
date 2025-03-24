## Superparamagnetism

`````{admonition} Source of text
:class: note
This text on superparamagnetic behavior at low temperature is from the Fall 2009 issue of the IRM Quarterly: *Bowles, Julie; Jackson, Mike; Chen, Amy; Solheid, Peat. (2009). IRM Quarterly, Volume 19, Number 3 (Fall 2009). Cover article: Interpretation of Low-Temperature Data Part 1: Superparamagnetism and Paramagnetism. Retrieved from the University Digital Conservancy, https://hdl.handle.net/11299/171316.*
 with minor edits by Nick Swanson-Hysell.
`````

Superparamagnetism (SP) describes the state of a single-domain-sized grain when thermal energy is sufficient to overcome barriers to a reversal of magnetization. Here the term “reversal” assumes uniaxial anisotropy with two minimum energy states having anti-parallel moment orientations. Energy barriers to magnetization changing arise from magnetocrystalline, magnetoelastic and/or shape anisotropy, all of which are proportional to grain volume (V). When the energy barriers are large with respect to thermal energy, the magnetization is “blocked” and the probability of spontaneous reversal approaches nil. But when the barriers are relatively low, thermal excitations can result in reversal of the magnetization over very short time scales, and the grain is in a superparamagnetic state. At a given temperature the volume at which a particle goes from being unblocked to blocked is known as the blocking volume (V<sub>b</sub>). For a given volume, we can block the grain by lowering the temperature (i.e. decreasing the available thermal energy) below the blocking temperature (T<sub>b</sub>).

```{figure}
https://upload.wikimedia.org/wikipedia/commons/thumb/d/d8/Ferrofluid_close.jpg/800px-Ferrofluid_close.jpg
:name: ferrofluid_figure
:align: center
:width: 60%

Ferrofluid (containing SP iron particles) placed over a
rare earth magnet. Source: [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Ferrofluid_close.jpg).
```

More formally, we can think about SP blocking in terms of a particle’s characteristic relaxation time (τ) and how rapidly a particle or assemblage of particles may approach equilibrium. Based on Néel theory of thermally activated magnetization (Néel, 1949) and the terminology of Dunlop and Özdemir (1997), the equilibrium magnetization \(M\) for aligned, uniaxial grains is given by:

\[
M_{\mathrm{eq}}(H_{0}, T) \;=\; M_{s}\,\tanh(\alpha),
\]
where \(\alpha = \frac{\mu_{0}\,M_{s}\,H_{0}\,V}{k_{B}\,T}\). For small \(\alpha\),  
\[
M_{\mathrm{eq}} \;\approx\; M_{s}\,\alpha.
\]

The relaxation time \(\tau\) in weak or zero field is:
 
\[
\tau \;=\; \tau_{0}\,\exp\!\Bigl(\frac{K\,V}{k_{B}\,T}\Bigr),
\]
where \(M_{s}\) is the saturation magnetization, \(\mu_{0}\) is the permeability of free space, \(H_{0}\) is a (weak) applied field, and \(\tau_{0}\) is a time constant related to the atomic reorganization time (on the order of \(10^{-9}\) s). \(K\) is an anisotropy constant also related to particle microcoercivity. In the equilibrium state, some fraction of the population has the high‐energy moment orientation (antiparallel to the applied field) while others occupy the low‐energy orientation, and spontaneous reversals occur at the same rate in each orientation (like a two‐well system in chemical equilibrium). For an assemblage of randomly oriented single‐domain grains, the relationship follows a Langevin function as:

\[
M_{\mathrm{eq}}(H_{0}, T) \;=\; M_{s}\,\Bigl[\coth(\alpha) \;-\; \frac{1}{\alpha}\Bigr].
\] 

\[
\approx \frac{M_s \alpha}{3}
\]
for small \(\alpha\).

When \(\tau\) is short compared to the observation time (the superparamagnetic state), the sample moment can equilibrate quickly with changes in applied field. In contrast, when \(\tau\) is long the sample is blocked and in a stable single domain (SSD) state; it can take a very long time (millions or billions of years) to reach magnetic equilibrium. Whether or not a particle is SP or SSD will depend on both temperature and volume, through their affect on \(\tau\).

## Remanence Data

The typical remanence experiments performed by visitors to the IRM include low-temperature cycling (300K → 20K → 300K) of a room-temperature (300K) saturation isothermal remanence (SIRM). A second common experiment cools the sample in zero field (ZFC) to low (e.g. 20K) temperature, where an SIRM is acquired. This low-temperature remanence is then measured as the sample warms back to room temperature. A room-temperature SIRM is not typically useful for studying SP particles, because they will not block at room temperature. However, many grains that are SP at room temperature will be SSD (blocked) at 20K. A sample with a **low-temperature SIRM that warms in zero field** (\(M_{eq} = 0\)) will demagnetize via the unblocking of these particles as they pass through their blocking temperature. If the sample contains a relatively narrow grain-size distribution of SP particles, the low-T remanence will unblock over a relatively narrow temperature interval (e.g. Fig. 2, circles). If, however, the sample has a more distributed grain-size distribution, the unblocking will be more gradual (e.g. Fig. 2, squares and diamonds). See Worm and Jackson (1999) for a discussion of volume distribution calculations based on the unblocking of a low-temperature remanence.

```{figure}
https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/SP_ZFC.png
:name: SP_ZFC
:align: center
:width: 60%

Saturation IRM acquired at 10K (ZFC) and measured on warming from 10K to 300K. The decrease in magnetization results from thermal unblocking of an SP population. Blue circles are from a Tiva Canyon Tuff sample with a relatively narrow grain-size distribution. Data from the same sample are shown in Fig. 3 and Fig. 4b. Red squares and green triangles are from synthetic, glassy basalts with more distributed grain-size distributions.
```

Not every such gradual decrease of low-temperature temperature remanence on warming is related to superparamagnetism, however. Moskowitz et al. (1998) demonstrate that multi-domain, high-Ti titanomagnetite exhibits similar behavior resulting from domain reorganization as magnetocrystalline anisotropy decreases on warming. (This phenomenon will be discussed at greater length in a future *IRM Quarterly* article.) It has also been demonstrated that partially oxidized magnetite (for example) may exhibit similar behavior. (See, e.g., Visiting Fellow report by Amy Chen, pg. 4, Fig. 2, insets.)

## Induced Magnetization Data

In an experiment more common to the physics community (but increasingly used by rock magnetists, e.g. Berquó et al., 2007), a sample is cooled in zero field (**ZFC**) to 20 K (for example), and the induced sample moment (\(M_i\)) is then measured on warming in the presence of a small (e.g. 5 mT) field. **This is followed by cooling the sample in the same 5 mT field (FC) and again measuring \(M_i\) on warming in a 5 mT field (Fig. 3).** In the FC case, for \(T < T_b\), magnetization can be described by regular TRM theory, i.e. \[
M_{TRM}(T) = M_{eq}(T_B) \frac{M_S(T)}{M_S(T_B)}
\]  
**Eq. 4**  

for a uniform distribution of particles with a single \( T_B \). For non-uniform distributions (e.g. Fig. 3), Eq. 4 must be integrated over all \( T_B \).  
**In the ZFC case, as the sample warms back up,** \( M_{i(ZFC)} \) **is initially less than** \( M_{i(FC)} \) **because the sample has blocked in a zero-field remanence.** As the sample warms through \( T_B \), \( \tau \) becomes short enough that the moment approaches \( M_{eq} \) in the applied field. Where the FC and ZFC curves join together defines the (maximum) blocking temperature, and for \( T > T_B \), both curves represent \( M_{eq} \). If the sample has a narrow-enough grain-size distribution and enough is known about its properties in order to estimate \( K \), particle volume can be estimated by solving Eq. 3 for \( V \).  

## AC Susceptibility Data  

In measuring low-field susceptibility (\( \chi \)), a weak, time-varying alternating field \([H(t) = H_0 \cos(\omega t)]\) is applied to the sample, and the magnetic response is measured. Consider the behavior of a superparamagnetic grain above and below its blocking temperature.  

At \( T_B << T < T_C \), **the grain is ferromagnetically ordered**, but \( \tau \) is very short compared to the observation time (\( t_{obs} \)), and the moment will respond **exactly in phase** with the alternating field.  

For an assemblage of randomly-oriented identical grains, susceptibility falls off with \( 1/T \) according to:  

\[
\chi_{sp} = \frac{\mu_0 V M_S^2}{3kT} = \frac{M_{eq}(H_0)}{H_0}
\]  
**Eq. 5**  

This type of behavior can be observed in Fig. 4a at \( T > \sim200K \).  

At \( T << T_B \), \( \tau \) **is very long compared to** \( t_{obs} \), the grains are fully blocked in the SSD state, and susceptibility will be **negligible in weak fields**:  

\[
\chi_{sd} = \frac{2M_S}{3H_k}
\]  
**Eq. 6**  

where \( H_k \) is the microcoercivity of the grain (the field required to reverse the magnetization in the absence of any thermal effects). This type of behavior is observed in Fig. 4c at \( T < \sim100K \).





