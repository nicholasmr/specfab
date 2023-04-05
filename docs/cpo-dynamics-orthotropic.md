# CPO dynamics for orthotropic grains


## Introduction

This tutorial focuses on modelling the CPO evolution of polycrystalline olivine, i.e. distributions of slip-system axes $n(\theta,\phi)$ and $b(\theta,\phi)$.

The tutorial shows how to model the CPO evolution of a Lagrangian material parcel subject to three different modes of deformation/stress:

![](https://raw.githubusercontent.com/nicholasmr/specfab/main/images/deformation-modes/deformation-modes.png#center){: style="width:620px"}

### Notation

Given the expansions

$$
n({\bf x},t,\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}n_{l}^{m}({\bf x},t) Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip-plane normals)}, 
\\
b({\bf x},t,\theta,\phi)=\sum_{l=0}^{L}\sum_{m=-l}^{l}b_{l}^{m}({\bf x},t) Y_{l}^{m}(\theta,\phi) \quad\text{(distribution of slip directions)}, 
$$

CPO evolution can be written as a matrix problem involving

$$
{\bf s}_n = [n_0^0,n_2^{-2},n_2^{-1},n_2^{0},n_2^{1},n_2^{2},n_4^{-4},\cdots,n_4^{4},\cdots,n_L^{-L},\cdots,n_L^{L}]^{\mathrm{T}} \quad\text{($n$ state vector)},
\\
{\bf s}_b = [b_0^0,b_2^{-2},b_2^{-1},b_2^{0},b_2^{1},b_2^{2},b_4^{-4},\cdots,b_4^{4},\cdots,b_L^{-L},\cdots,b_L^{L}]^{\mathrm{T}} \quad\text{($b$ state vector)},
$$

such that 

$$
\frac{\mathrm{D}{\bf s}_n}{\mathrm{D} t} = {\bf M}_n \cdot {\bf s}_n \quad\text{($n$ state evolution)},
\\
\frac{\mathrm{D}{\bf s}_b}{\mathrm{D} t} = {\bf M}_b \cdot {\bf s}_b \quad\text{($b$ state evolution)},
$$

where the operators (matrix) ${\bf M}_n$ and ${\bf M}_b$ represents the effect of a given CPO process, which may depend on stress, strain-rate, temperature, etc.

!!! note
    Only lattice rotation is so far supported: ${\bf M} = {\bf M}_{\mathrm{LROT}}$ for both $n$ and $b$.

- - -

## Lattice rotation

To be published before documented here.

<!--
### Example
-->

- - -

## Regularization

Same as [CPO dynamics for transversely isotropic grains](cpo-dynamics-tranisotropic.md)



