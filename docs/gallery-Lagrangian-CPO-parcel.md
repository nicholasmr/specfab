# Lagrangian CPO parcel

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/modes-strain/lagrangian-parcel-trajectory.png#center){: style="width:400px"} 

A Lagrangian parcel refers to a small, moving material volume that is followed through time as it flows within e.g. a glacier or ice sheet, unlike the [Eulerian perspective](gallery-Eulerian-CPO-field.md) which focuses on fixed locations in space. 
A Lagrangian description is well-suited for studying CPO evolution along a flow line ${\bf x}(t)$ if the thermomechanical background conditions are (approximately) steady; that is, the velocity, temperature, and stress fields are constant in time, *though not necessarily in space*.
By treating each parcel as a discrete entity with its own microstructural state that updates with time, a Lagrangian parcel model is a particularly simply way to estimate CPO evolution and to calculate CPO-induced quantities along a flow line, such as mechanical or dielectric anisotropy. 

Below, some examples are given on how to model CPO evolution of a Lagrangian parcel relevant to glacier ice.

## Constant conditions

*Specfab* includes a high-level integrator for calculating the CPO evolution of a Lagrangian parcel subject to a spatio-temporally constant strain-rate, stress and temperature.
The following code illustrates how to use it, which relies on specifying the kinematic mode of deformation in terms of the `DK` object.

```python
--8<-- "docs/snippets/Lagrangian-CPO-parcel/constant-conditions.py"
```

## Ice core CPOs

A Lagrangian approach is well-suited for modelling the vertical CPO profile at ice sheet domes and divides. 
Assuming, for example, the classical Nye model of an ice divide of height $H$ (no basal melt, constant rate of thinning, a constant accumulation rate $a$), the [velocity gradient](deformation-kinematics.md) is constant and equal to 

$$
\nabla {\bf u} = 
-\frac{1}{\tau}
\begin{bmatrix}
-(1+q)/2 & 0 & 0\\
0 & -(1-q)/2 & 0\\
0& 0& 1
\end{bmatrix}
,
$$

where $q=0$ for a axis-symmetric dome and $q=\pm 1$ for a divide aligned with the $x$ or $y$ direction, and the $e$-folding time $\tau$ is 

$$    
\tau = \frac{a}{H}.
$$

###📝 Code example
The example shows how to model the CPO profile of the GRIP ice core, Greenland:

```python
--8<-- "docs/snippets/Lagrangian-CPO-parcel/GRIP.py"
```

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/snippets/Lagrangian-CPO-parcel/GRIP.py#center){: style="width:300px"} 

## SSA parcel

**To be documented...**
