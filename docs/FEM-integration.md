# Finite element integration

## The coupled problem 

The coupled evolution of flow and CPO can be solved by joining different components of specfab together with a Stokes flow solver (purple piece):

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/jigsaw/jigsaw1.png#center){: style="width:360px"}

That is, specfab can model/provides

* [CPO evolution](cpo-dynamics-tranisotropic.md) given the local large-scale stress, strain-rate and temperature (green piece)
* [Viscous anisotropies](enhancements-strainrate.md) induced by the local CPO (orange piece)
* [Bulk anisotropic rheologies](constitutive-viscoplastic.md) given the local viscous anisotropy (blue piece)

## Elmer/Ice

See [Lilien et al. (2023)](https://doi.org/10.1017/jog.2023.78) for ice-flow modelling.

## FEniCS

To be updated.	


<!--
## Icepack

To be updated.	

## Úa

To be updated.	

-->


