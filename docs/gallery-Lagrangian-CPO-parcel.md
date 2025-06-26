# Lagrangian CPO parcel

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/modes-strain/lagrangian-parcel-trajectory-ice.png#center){: style="width:380px"} 

A Lagrangian parcel refers to a small, moving material volume that is followed through time as it flows within e.g. a glacier or ice sheet, unlike the [Eulerian perspective](gallery-Eulerian-CPO-field.md) which focuses on fixed locations in space. 
A Lagrangian description is well-suited for studying CPO evolution along a flow line ${\bf x}(t)$ if the thermomechanical background conditions are approximately steady; that is, the velocity, temperature, and stress fields are constant in time, *though not necessarily in space*.
By treating each parcel as a discrete entity with its own microstructural state that updates with time, a Lagrangian parcel model is a particularly simply way to estimate CPO evolution and to calculate CPO-induced quantities along a flow line (trajectory), such as mechanical or dielectric anisotropy. 

Below, some examples are given on how to model CPO evolution using the Lagrangian parcel approach, relevant to glacier ice.

- - -

## Ice divide

/// html | div[style='float: left; width: 63%;']
A Lagrangian approach is well-suited for modelling the vertical CPO profile at domes and divides of ice sheets. 
Assuming e.g. the classical Nye model of an ice divide of height $H$ (no basal melt, constant rate of thinning, a constant accumulation rate $a$), the [velocity gradient](deformation-kinematics.md) is constant and equal to 

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

where $q=0$ for a axis-symmetric dome and $q=\pm 1$ for a divide aligned with the $x$ or $y$ direction. 
The $e$-folding time $\tau$ is a function of the ice-equivalent accumulation rate and divide thickness: 

$$    
\tau = \frac{a}{H}.
$$

Depending on whether temperatures are large enough to [activate DDRX](fabdyn-DDRX.md) or not, the temperature profile must be prescribed, too. 

///

/// html | div[style='float: right;width: 2%;']
///

/// html | div[style='float: right;width: 35%;']
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation/divide-parcel.png){: style="width:300px"} 
///

/// html | div[style='clear: both;']
///

###📝 Example: GRIP ice core

/// html | div[style='float: left; width: 38%;']
The following code example shows how to model the CPO profile of the GRIP ice core, Greenland. 
Eigenvalues of other ice cores can be [found here](https://github.com/nicholasmr/specfab/tree/main/data/icecores).

The figure to the right shows the model result (lines) compared to observations (markers) from thin-section analysis. 
At around $z/H \simeq 0.2$, temperatures are sufficiently high that DDRX becomes important and affects (weakens) the vertical single-maximum CPO that otherwise strengthens with depth.  
///

/// html | div[style='float: right;width: 2%;']
///

/// html | div[style='float: right;width: 60%;']
![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/snippets/Lagrangian-CPO-parcel/GRIP-parcel.png#center){: style="width:400px"} 
///

/// html | div[style='clear: both;']
///

```python
--8<-- "docs/snippets/Lagrangian-CPO-parcel/GRIP-parcel.py"
```

- - -

## Plug flow

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/SSA-fabric/lagrangian-column-trajectory.png#center){: style="width:450px"} 

If the Shallow Shelf/Stream Approximation (SSA) is applicable, velocities can be assumed depth constant (i.e., *plug flow* without vertical shearing). 
In this case, a [depth-average treatment of CPO evolution](gallery-Eulerian-CPO-field.md) transforms the Lagrangian parcel model into a Lagrangian *column* model. 
This generalizes the above parcel model, since the velocity gradient and stress fields can no longer be assumed constant but depend on the column position ${\bf x}(t)=[x(t),y(t)]$. 

Assuming that both velocity and CPO fields are steady, and neglecting surface accumulation that adds new (isotropic) ice, the problem can be solved by prescribing satellite-derived surface velocities and estimates of the ice temperature.  

###📝 Example: Pine Island Glacier

The following code example shows how to model CPO evolution along a flow line over the Pine Island Glacier, Antarctica. 
The example relies on *MEaSUREs* ice velocities, assumes a constant temperature field, and asserts that the stress tensor can be approximated as coaxial to the strain rate tensor to close the DDRX problem. 
If the ice mass is sufficiently cold, DDRX can be neglected altogether and only satellite-derived velocities are needed to close the problem. 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/snippets/Lagrangian-CPO-parcel/PIG-column.png#center){: style="width:700px"} 

```python
--8<-- "docs/snippets/Lagrangian-CPO-parcel/PIG-column.py"
```

