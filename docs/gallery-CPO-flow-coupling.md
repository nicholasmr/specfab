# Coupling flow and CPO

The coupled evolution of flow and CPO can be solved by joining different components of specfab together with a momentum balance solver (e.g., Stokes solver):

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/jigsaw/jigsaw2.png#center){: style="width:340px"}

That is, specfab can model:

* [CPO evolution](fabdyn-matrix-model.md) given bulk stress, strain-rate and temperature fields (green piece)
* [Viscous anisotropy](enhancements-strainrate.md) induced by the CPO (orange piece)
* [Bulk anisotropic rheologies](constitutive-viscoplastic.md) given the local viscous anisotropy (blue piece)

### Elmer/Ice

See [Lilien et al. (2023)](https://doi.org/10.1017/jog.2023.78) for ice flow modelling.

### FEniCS

See [Rathmann et al. (2024)](https://doi.org/10.1029/2024GC011831) for olivine modelling.

<!--
## Icepack

To be updated.	

## Úa

To be updated.	

-->


