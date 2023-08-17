# Deformation modes

For a continuum subject to deformation, the [deformation gradient tensor](https://www.continuummechanics.org/deformationgradient.html), ${\bf F}$, describes the relative change in position of material points.
If ${\bf F}$ is known, then the [velocity gradient tensor](https://www.continuummechanics.org/velocitygradient.html) follows as
$$
\nabla {\bf u} = \dot{{\bf F}} {\bf F}^{-1}
.
$$

The resulting strain experienced (strain tensor) is
\begin{align}
{\boldsymbol \epsilon} = \frac{1}{2}\left( {\bf F}+{\bf F}^\top \right) - {\bf I}
.
\end{align}

Below, we consider how to represent pure shear and simple shear with ${\bf F}$.

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/deformation-modes-2.png#center){: style="width:690px"}

## Pure shear

${\bf F}$ is diagonal for pure shear deformation when the principal strain axes are aligned with the coordinate system. 
Suppose shortening takes places along the vertical axis and lengthening in the horizontal plane, then 

$$
{\bf F}_{\mathrm{P}} = 
\begin{bmatrix}
b^{(1+r)/2} & 0 & 0\\
0 & b^{(1-r)/2} & 0\\
0& 0& b^{-1}
\end{bmatrix}
,
$$

where the parameter $r\in[-1;1]$ controls the relative lengthening between two horizontal directions: 

* for $r=0$ lengthening is equal in the $x$ and $y$ directions,

* for $r=+1$ lengthening occurs only in the $x$ direction,

* for $r=-1$ lengthening occurs only in the $y$ direction.

Calculating the velocity gradient tensor yields

$$
\nabla {\bf u} = 
\frac{\dot{b}}{b}
\begin{bmatrix}
(1+r)/2 & 0 & 0\\
0 & (1-r)/2 & 0\\
0& 0& -1
\end{bmatrix}
.
$$

If the scaling parameter, $b$, is written in terms of the $e$-folding time scale $T$ as 

$$
b(t) = \exp(t/T)
,
$$

the strain-rate becomes constant,

$$
\frac{\dot{b}}{b} = \frac{1}{T}. 
$$

Notice that the vertical strain experienced as a function of time is 

$$
\epsilon_{zz}(t) = \frac{1}{b(t)} - 1,
$$

where $\epsilon_{zz} = 0$ corresponds to an undeformed ice parcel ($t=0$), 
and $\epsilon_{zz} = -1$ is the limit of vanishing parcel height ($t\rightarrow\infty$). 

Consider instead the case where lengthening takes place along the vertical axis and shortening in the horizontal plane. 
This is achieved by substituting $b\rightarrow b^{-1}$ which implies ${\bf F}={\bf F}_{\mathrm{P}}^{-1}$ and hence is the time reversed behavior of ${\bf F}_{\mathrm{P}}$ since $b^{-1}(t)=b(-t)$.

### Example 

The above expressions are accessible in specfab as follows

```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

axis  = 2 # axis of shortening (T>0) or lengthening (T<0): 0=x, 1=y, 2=z
Tc    = 1 # time taken in seconds for parcel to reduce to half (50%) height if T>0, or abs(time) taken for parcel to double in height (200%) if T<0.
r     = 0 # asymmetry parameter for shortening (if T>0) or lengthening (if T<0)

T = Tc/np.log(2) # corresponding e-folding time
ugrad = sf.pureshear_ugrad(axis, r, T) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)      # strain-rate and spin tensor

t = 1 # some specific time of interest
b   = sf.pureshear_b(T, t)          # scaling parameter at time t
F   = sf.pureshear_F(axis, r, T, t) # deformation gradient tensor at time t
eps = sf.F_to_strain(F)             # strain tensor at time t
```

## Simple shear

Simple shear strain may be characterized by the shear angle $\gamma$ of the resulting rhombus.
In the case of vertical shear in the $x$&mdash;$z$ plane, the deformation gradient tensor is given by

$$
{\bf F}_{\mathrm{S}}
=
\begin{bmatrix}
1 & 0 & \tan(\gamma) \\
0 & 1 & 0\\
0& 0& 1
\end{bmatrix}
,
$$

and velocity-gradient tensor becomes

$$
\nabla {\bf u} = 
\frac{\dot{\gamma}}{\cos^2(\gamma)}
\begin{bmatrix}
0 & 0 & 1\\
0 & 0 & 0\\
0& 0& 0
\end{bmatrix}
.
$$

For a constant shear rate, $1/T$, the shear-angle time dependence is
\begin{align}
\gamma(t) = \tan^{-1}(t/T)
,
\end{align}
where $T$ is the characteristic time taken to reach a shear of 1 from an undeformed state.
In this case, $\nabla {\bf u}$ is constant, too, since 

$$
\frac{\dot{\gamma}}{\cos^2(\gamma)} = \frac{1}{T}
.
$$

### Example 

The above expressions are accessible in specfab as follows

```python
import numpy as np
from specfabpy import specfab as sf
lm, nlm_len = sf.init(4) 

plane = 1 # plane of shear: 0=yz, 1=xz, 2=xy
T     = 1 # time taken in seconds for parcel to a reach shear strain of 1 (45 deg shear angle)

ugrad = sf.simpleshear_ugrad(plane, T) # velocity gradient
D, W = sf.ugrad_to_D_and_W(ugrad)      # strain-rate and spin tensor

t = 1 # some specific time of interest
gamma = sf.simpleshear_gamma(T, t)    # shear angle at time t
F     = sf.simpleshear_F(plane, T, t) # deformation gradient tensor at time t
eps   = sf.F_to_strain(F)             # strain tensor at time t
```

