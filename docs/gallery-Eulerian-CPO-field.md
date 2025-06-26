# Eulerian CPO field

If the velocity field ${\bf u}({\bf{x}},t)$ and CPO field ${\bf s}({\bf{x}},t)$ can be assumed steady, CPO evolution reduces to a high-dimensional boundary value problem in ${\bf s}({\bf{x}},t)$ that can easily be solved using e.g. the finite element method. 

When considering the flow glacier ice, it is common to simplify the problem by invoking the Shallow Shelf Approximation (SSA), a depth-integrated version of the full Stokes equations. 
This approximation conveniently reduces the problem to a two-dimensional horizontal, membrane-like flow with negligible vertical shear. 

Performing a similar calculation, the depth-average expression for steady CPO evolution takes the form [(Rathmann et al., 2025)](https://eartharxiv.org/repository/view/8861/) 

$$
({\bf u} \cdot \nabla) {\bf \bar s} = 
({\bf M}_{\mathrm{LROT}}+{\bf M}_{\mathrm{DDRX}}+{\bf M}_{\mathrm{CDRX}}) \cdot {\bf \bar s} 
+ \frac{a_{\mathrm{sfc}}}{H}( {\bf s}_{\mathrm{sfc}} - {\bf \bar s} ) 
+ \frac{a_{\mathrm{sub}}}{H}( {\bf s}_{\mathrm{sub}} - {\bf \bar s} ) 
,
$$

where ${\bf \bar s}(x,y)$ is the depth-average CPO state vector field, ${\bf{u}}(x,y)=[u_x(x,y),u_y(x,y)]$ is the horizontal surface velocity field, and $H$ is the ice thickness. 
The first term represents CPO advection along stream lines, and the second term represents the depth-average effect of crystal processes. 
The third and fourth terms are state-space attractors, causing ${\bf \bar s}$ to tend towards the characteristic CPO states of ice that accumulates on the surface ${\bf s}_{\mathrm{sfc}}$ or subglacially ${\bf s}_{\mathrm{sub}}$, depending on the positively-defined ice-equivalent accumulation rates $a_{\mathrm{sfc}}$ and $a_{\mathrm{sub}}$. 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/SSA-fabric/SSA-fabric-long.png){: style="width:700px"}

## Closure

If CPO development is dominated by lattice rotation (typical for cold ice), the problem is closed by specifying the horizontal surface velocity field (e.g., satellite-derived velocities), together with accumulation rates and the characteristic CPO state of accumulated ice (typically isotropic). 

If DDRX is non-negligible (typical for warm ice), the temperature and stress field must additionally be prescribed. 
In this case, the problem can be closed by assuming some temperature field and assert that the stress and strain-rate tensors are coaxial so that ${\boldsymbol\tau} \propto \dot{\boldsymbol\epsilon}$, although this has some limitations [(Rathmann et al., 2025)](https://eartharxiv.org/repository/view/8861/).


## Regularization 

Noise in surface velocity products, in addition to uncertainties due to model assumptions, can render the steady SSA CPO problem ill-posed. 
To solve this, Laplacian regularization of the form $\nu \nabla^2 {\bf \bar s}$ can be adding to the right-hand side of the problem, at the expense of limiting how large spatial CPO gradients are permitted. 
The strength of regularization $\nu$ (`nu_real` in code below) is therefore a free model parameter which must be carefully selected, especially in very dynamic regions where the CPO field might change rapidly with distance. 

## 📝 Example: Pine Island Glacier

Let us consider the Pine Island Glacier (PIG), Antarctica, as an example for how to solve the steady SSA CPO problem, assuming that: 

1. SSA is appropriate for both floating and grounded ice over the region.

2. CPO development is dominated by englacial crystal processes so that contributions from accumulation can be neglected. 

3. Ice temperatures can be approximated as constant over the domain (needed only if modeling DDRX). 

➡️ To begin, the model domain needs to be meshed. For this, we use *gmsh* with the `mesh.geo` file: 

```gmsh
--8<-- "docs/snippets/steady-SSA-solver/mesh.geo"
```
which looks like this:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/snippets/steady-SSA-solver/mesh.png){: style="width:480px"}

➡️ Next, the CPO problem is solved using *FEniCS*+*specfab* and results plotted as follows:

```python
--8<-- "docs/snippets/steady-SSA-solver/PIG.py"
```
The resulting plots are shown below, where the first row shows the velocity and strain rate field interpolated onto the model mesh, and the second and third row shows model results without and with DDRX, respectively. 

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/docs/snippets/steady-SSA-solver/PIG-gallery.png){: style="width:700px"}


