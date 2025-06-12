# Deformation kinematics

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation/displacement-field.png#center){: style="width:630px"}

The [deformation gradient tensor](https://www.continuummechanics.org/deformationgradient.html) ${\bf F}$ quantifies how infinitesimal material elements (line, area, or volume) are deformed relative to some original reference configuration.
The strain tensor ${\boldsymbol \epsilon}$ is derived from ${\bf F}$ and is a measure of how a material deforms, excluding rigid body motions like rotation and translation.
In continuum mechanics, it is useful to isolate such pure deformation (stretching) from rigid body motion, since the viscous and elastic response of a material when subject external forces is concerned with how easy or hard it is to stretch a material.

## Small strain tensor 

The symmetric and antisymmetric parts of ${\bf F}$ describe the degree of stretching and the rotation of material elements, respectively, useful for isolating pure deformation (stretching) from rigid body rotations. 
For small deformations, the strain tensor (formally *small* strain tensor) is therefore defined as the symmetric part of ${\bf F}$, minus the identity:
$$
{\boldsymbol \epsilon} = \frac{1}{2}\left( {\bf F}+{\bf F}^\top \right) - {\bf I}.
$$
The strain tensor ${\boldsymbol \epsilon}$ characterizes how a material stretches, compresses, or shears: positive components indicate extension along certain directions, negative components indicate compression, and off-diagonal components correspond to shear deformation.

The antisymmetric part of ${\bf F}$ is related to the component of deformation that describes rigid body rotation (so-called rotation tensor), defined as
$$
{\bf w} = \frac{1}{2}\left( {\bf F}-{\bf F}^\top \right) + {\bf I}.
$$
Notice that ${\bf F}= {\boldsymbol \epsilon} + {\bf w}$ by definition.

While the strain tensor is the kinematic field of interest for elastic deformation (elastic problems are concerned with displacement), the strain-rate and spin tensors are more relevant for viscous deformation (viscous problems are concerned with velocity, ${\bf u}$). 
The strain-rate and spin tensors are simply the rate-of-change of ${\boldsymbol \epsilon}$ and ${\bf w}$, respectively:
$$
{\dot{\boldsymbol\epsilon}} = \frac{1}{2}\left(\nabla {\bf u} + \nabla {\bf u}^\top\right) ,
\quad
{\boldsymbol\omega} = \frac{1}{2}\left( \nabla {\bf u} - \nabla {\bf u}^\top \right) ,
$$
where the [velocity gradient tensor](https://www.continuummechanics.org/velocitygradient.html) is $\nabla {\bf u} = \dot{{\bf F}}  \cdot {\bf F}^{-1}$.

## Kinematic modes

Let us now turn to the two common kinematic modes of deformation, pure shear and simple shear, that we will encounter many times throughout this book.
Introducing these kinematic modes in terms of their deformation gradient ${\bf F}$ provides a powerful way to describe e.g. how ice deforms below ice divides and domes of ice sheets.

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/modes-strain/kinematic-modes.png#center){: style="width:620px"}


!!! warning "Example of strain tensors"
    
    ![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/modes-strain/strain-tensors.png#center){: style="width:570px"}


### Pure shear

Pure shear (stretching) aligned with the Cartesian coordinate axes renders ${\bf F}$ diagonal, a convenient frame to consider. 
Suppose shortening takes places along the vertical $z$-axis and lengthening in the horizontal $xy$ plane so that volume is conserved ($\det({\bf F})=1$; incompressible material). The deformation gradient can then be written as 

$$
{\bf F} = 
\begin{bmatrix}
r^{-(1+q)/2} & 0 & 0\\
0 & r^{-(1-q)/2} & 0\\
0& 0& r
\end{bmatrix}
,
$$

where $r$ is a scaling parameter that controls the relative vertical shortening of material lines, and $q\in[-1;1]$ determines the horizontal direction of lengthening:

* For $q=0$ lengthening is equal in the $x$ and $y$ directions. 
* For $q=+1$ lengthening occurs only in the $x$ direction. 
* For $q=-1$ lengthening occurs only in the $y$ direction.

The corresponding vertical strain with time follows immediately as 
\begin{align}
\epsilon_{zz}(t) = r(t) - 1,
\end{align}
where $\epsilon_{zz} = 0$ (or $F_{zz}=1$) corresponds to the undeformed configuration and $\epsilon_{zz} = -1$ (or $F_{zz}=0$) corresponds to the limit of a vanishing material (parcel) height.
Conversely, if lengthening takes place along the vertical axis with shortening in the horizontal plane then ${\bf F}\rightarrow {\bf F}^{-1}$ above, implying $r\rightarrow r^{-1}$.

**Constant strain rate**<br>
For pure shear of the above form, it can be shown that the velocity gradient tensor reduces to 

$$
\nabla {\bf u} = 
\frac{\dot{r}}{r}
\begin{bmatrix}
-(1+q)/2 & 0 & 0\\
0 & -(1-q)/2 & 0\\
0& 0& 1
\end{bmatrix}
.
$$

Clearly, if $r$ is written in terms of the $e$-folding time $\tau$, 

$$
r(t) = \exp(-t/\tau) ,
$$ 

the strain-rate tensor becomes constant. Specifically, the vertical component becomes

$$
\epsilon_{zz} = -\frac{1}{\tau} .
$$

**Code example**<br>
The above expressions are accessible in specfab as follows

```python
--8<-- "docs/snippets/defkin-pureshear.py"
```

### Simple shear

Simple shear may be characterized by the shear angle $\gamma$ of the resulting rhombus.
In the case of vertical shear in the $xz$ plane, the deformation gradient tensor is given by

$$
{\bf F} =
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

**Constant strain rate**<br>
For a constant rate of shear $1/\tau$, where $\tau$ is a characteristic time, the prefactor is  

\begin{align}
\frac{\dot{\gamma}}{\cos^2(\gamma)} = \frac{1}{\tau} .
\end{align}

Solving for the time-dependent shear angle function can be shown to give

$$
\gamma(t) = \tan^{-1}(t/\tau) .
$$

Thus, $\tau$ is the time it takes to reach a shear angle of $\gamma=45^\circ$ from the undeformed configuration, corresponding to $\epsilon_{xz}=\tan(45^\circ)/2=0.5$.

**Code example**<br>
The above expressions are accessible in specfab as follows

```python
--8<-- "docs/snippets/defkin-pureshear.py"
```

- - -

## Plotting deformed parcels

Given a deformation gradient ${\bf F}$, the resulting deformed parcel shape can be plotting as follows: 

```python
--8<-- "docs/snippets/defkin-parcel.py"
```
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/tests/deformation-kinematics/deformed-parcel.png){: style="width:450px"}

