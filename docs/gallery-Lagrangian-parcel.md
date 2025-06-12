# Lagrangian parcel model

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/modes-strain/lagrangian-parcel-trajectory.png#center){: style="width:400px"} 

*Specfab* includes a high-level integrator for calculating the CPO evolution of a Lagrangian parcel subject to a constant strain rate tensor.
The following code illustrates how to use it, which relies on specifying the kinematic mode of deformation in terms of the `DK` object.

```python
--8<-- "docs/snippets/gallery-lagrangian-parcel.py"
```

