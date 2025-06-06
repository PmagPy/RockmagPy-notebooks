# Introduction to anisotropy of magnetic susceptibility (AMS)

Magnetic susceptibility describes the relationship between an applied magnetic field and the magnetization that is induced by this field:

$M = \chi H$

where $M$ is the magnetization (in Am<sup>2</sup>/kg), $H$ is the applied field (in A/m), and $\chi$ is the magnetic susceptibility (this is the mass susceptibility in m<sup>3</sup>/kg). Magnetic susceptibility reflects the diamagnetic, paramagnetic, and/or ferromagnetic response of minerals within a sample.

If a material is perfectly isotropic with the magnetization response being uniform in all directions, the relationship described by $\chi$ is independent of the orientation of the sample relative to the applied field. In this case, $\chi$ is a scalar quantity (a single number). However, if there is a different response in different orientations the material can be considered to be anisotropic. In this case, $\chi$ needs to be described as a tensor ($\chi_{ij}$).

In practice, anisotropy of magnetic susceptibility (AMS) is quantified by varying the orientation of a sample within the coils of a magnetic susceptibility bridge. By doing so, the length and orientation of the principal, major, and minor eigenvectors which are defined as:

$V_1 \geq V_2 \geq V_3$

these eigenvectors are also referred to in the literature as:

$K_{max} \geq K_{int} \geq K_{min}$

and the AMS ellipsoid can be illustrated as:

<figure>
  <img src="https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/AMS.png" alt="A three-dimensional ellipsoid with orthogonal principal axes labeled Kmax, Kint, and Kmin, representing the maximum, intermediate, and minimum magnetic susceptibility directions. A plane intersects the Kmax and Kint axes and is labeled magnetic foliation." width="450">
  <figcaption>From https://www.agico.cz/downloads/documents/agicoprints/ams_rocks_pres.pdf </figcaption>
</figure>

In the above plot, the axes of the ellipsoid are the eigenvectors (commonly referred to as the “principal axes”) with their magnitudes corresponding to the eigenvalues. Collectively, the eigenvalues and eigenvectors are called the eigenparameters. The orientation of this ellipsoid and the relative magnitude of the eigenvectors is central to the interpretation of AMS data. 

## The MagIC convention for the magnetic susceptibility tensor

The `aniso_s` list is a six-element colon-delimited list provided in the `aniso_s` column of a MagIC specimens tables that is a succinct representation of the anisotropy tensor $\chi$. In the linear equations below, $M_1$, $M_2$, and $M_3$ are the components of the magnetization vector $M$. The relationship between the magnetization and the magnetic field ($H_1$, $H_2$, and $H_3$) is described by the components of the magnetic susceptibility tensor $\chi_{ij}$. If the material is anisotropic, there will be different values in different orientations described by the $\chi$ tensor.

$ M_1 = \chi_{11} H_1 + \chi_{12} H_2 + \chi_{13} H_3 $

$ M_2 = \chi_{21} H_1 + \chi_{22} H_2 + \chi_{23} H_3 $

$ M_3 = \chi_{31} H_1 + \chi_{32} H_2 + \chi_{33} H_3 $

This susceptibility tensor has six independent matrix elements because $\chi_{ij} = \chi_{ji}$. The `aniso_s` column contains these six elements as $s_1$ : $s_2$ : $s_3$ : $s_4$ : $s_5$ : $s_6$ which are defined as:

$ s_1 = \chi_{11} $

$ s_2 = \chi_{22} $

$ s_3 = \chi_{33} $

$ s_4 = \chi_{12} = \chi_{21} $

$ s_5 = \chi_{23} = \chi_{32} $

$ s_6 = \chi_{13} = \chi_{31} $

## Anisotropy measurements in practice

In many laboratories, magnetic susceptibility is measured by placing a sample in a solenoid with an applied field *H*. The induced magnetization *M* parallel to *H* is measured in different orientations. Only s1, s2, and s3 can be measured directly. The other terms s4, s5, and s6 are only indirectly determined. For multi-positional approaches, the sets of values of measurements such as susceptibility would be $K_i = A_{ij}s_j$ where *A* is the design matrix. 

### The Jelínek 15 position diagram

**Statistical processing of anisotropy of magnetic susceptibility measured on groups of specimens** by Jelínek (1978) illustrates the 15 position scheme of anisotropy measurement. 

<figure>
  <img src="https://raw.githubusercontent.com/PmagPy/RockmagPy-notebooks/main/book/images/Jelinek_15_method.png" alt="15 position scheme of anisotropy measurement" width="600">
</figure>

It is also possible to adopt a 9-position experimental approach. In this case, follow positions 1,2,3, 6,7,8, 11,12,13 of the Jelínek 15 position scheme for the anisotropy experiments. 

The 15-position design matrix is given by:
$
A = \begin{pmatrix}
.5 & .5 & 0 & -1 & 0 & 0 \\
.5 & .5 & 0 & 1  & 0 & 0 \\
1  & 0  & 0 & 0  & 0 & 0 \\
.5 & .5 & 0 & -1 & 0 & 0 \\
.5 & .5 & 0 & 1  & 0 & 0 \\
0 & .5 & .5 & 0  & -1 & 0 \\
0 & .5 & .5 & 0  & 1  & 0 \\
0 & 1 & 0 & 0  & 0 & 0 \\
0 & .5 & .5 & 0  & -1 & 0 \\
0 & .5 & .5 & 0  & 1  & 0 \\
.5 & 0 & .5 & 0  & 0 & -1 \\
.5 & 0 & .5 & 0  & 0 & 1 \\
0 & 0 & 1 & 0  & 0 & 0 \\
.5 & 0 & .5 & 0  & 0 & -1 \\
.5 & 0 & .5 & 0  & 0 & 1 \\
\end{pmatrix}
$

This 15-position procedure is widely adopted in experiments that seek to quantify the AMS tensor.