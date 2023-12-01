#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

import copy, sys, time, code # code.interact(local=locals())

import numpy as np
from scipy.spatial.transform import Rotation as R

from dolfin import *
from specfabpy.fenics.rheology import Orthotropic, Isotropic

### Mesh

N = 5
mesh = UnitCubeMesh(N,N,N)
plot(mesh, title="Unit cube")

T = FunctionSpace(mesh, FiniteElement('CG', mesh.ufl_cell(), 1))
U = FunctionSpace(mesh, VectorElement('CG', mesh.ufl_cell(), 1))
W = FunctionSpace(mesh, TensorElement('CG', mesh.ufl_cell(), 1))

### mi

r0 = R.from_euler('z', 68, degrees=True)
r1 = R.from_euler('x', 14, degrees=True)
r = r0*r1
rotmat = r.as_matrix()
m1,m2,m3 = np.eye(3) 
m1,m2,m3 = np.matmul(rotmat,m1), np.matmul(rotmat,m2), np.matmul(rotmat,m3)
mi  = [Constant(m1), Constant(m2), Constant(m3)]

#print(m1, m2, m3, np.linalg.norm(m1), np.linalg.norm(m2), np.linalg.norm(m3))

### Eij

A   = Constant(1)

#Eij = [Constant(1)]*6
Eij = [Constant(0.95), Constant(0.7), Constant(0.8),  Constant(2.5), Constant(4), Constant(10)] # random enhancements

### Strain-rate tensors

u = Expression(('u0*x[2]','0','0'), u0=100, degree=1)
D_ss_xz = sym(grad(project(u, U)))

u_ = Expression(('0','0','u0*x[1]'), u0=100, degree=1)
D_ss_yz = sym(grad(project(u, U)))

u_ = Expression(('0','u0*x[0]','0'), u0=100, degree=1)
D_ss_xy = sym(grad(project(u, U)))

D0 = 2.0
D_uc_xx = as_tensor(D0*np.diag([-2, +1, +1]))
D_uc_yy = as_tensor(D0*np.diag([-2, +1, +1]))
D_uc_zz = as_tensor(D0*np.diag([+1, +1, -2]))

### Loop over cases

n_list = [1,3,4]
D_list = [D_ss_xz, D_ss_yz, D_ss_xy, D_uc_xx, D_uc_yy, D_uc_zz, ]

D_list_names = ['SS xz', 'SS yz', 'SS xy', 'UC xx', 'UC yy', 'UC zz']

def error(D_test, D_true):
    err   = project(D_test-D_true, W)
    errsq = project(inner(err,err), T)
    errsqtot = assemble(errsq*dx)
    return errsqtot
    
for n in n_list:

    orth = Orthotropic(n=n)
    iso  = Isotropic(n=n)

    for ii, D_true in enumerate(D_list):

        # Test orthotropic forward--inverse pair        
        S_test = orth.S(D_true, A, mi, Eij)
        D_test = orth.D(S_test, A, mi, Eij)
        errsqtot = error(D_test, D_true)

        # Isotropic reference
        Siso_test = iso.S(D_true, A)
        Diso_test = iso.D(Siso_test, A)        
        errsqtot_iso = error(Diso_test, D_true)
        
        print('n=%i, D_true=%s :: sum(abs(error_ortho)) = %.1e  (iso err = %.1e)'%(n, D_list_names[ii], errsqtot, errsqtot_iso)) # should be close to zero

