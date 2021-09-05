#!/usr/bin/python3
# Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021

# Test isotropic and anisotropic rheologies as a free-standing cube of ice is subject to a compressional or shear traction force on its top face.
# Anisotropic parameters are set manually (see "Fabric" section).

import copy, sys, code # code.interact(local=locals())
import numpy as np
from dolfin import *
import matplotlib.pyplot as plt
from vedo.dolfin import plot as vplot

x1 = Constant((1,0,0))
x2 = Constant((0,1,0))
x3 = Constant((0,0,1))

# Constants
nglen  = 1 # Flow exponent
Aglen  = 1 # Rate factor (=1 is fine for these kind of experiments)
rhoi   = 917
gmag   = 0 # Deformation test without gravity
g      = Constant((0,0,-gmag)) # grav vector  

#------------------
# Deformation type
#------------------

#dtype = 'comp'
dtype = 'shear'

taushearmax = 1e1
if dtype=='shear': r =  x1 # "r" is the traction vector direction on top surface
if dtype=='comp':  r = -x3

y0_plane = 0.5 # y-coordinate of x--z plane to plot velocities components on

#------------------
# Fabric
#------------------

# Orthotropic directional enhancement factors

m1,m2,m3 = x1,x2,x3 # Symmetry axes (principal fabric axes)

E11 = 1 
E22 = 1
E33 = 1 # z compressional enhancement factor
E31 = 2 # x--z shear enhancement factor
E23 = 1
E12 = 1

Eij = Constant(np.matrix(((E11,E12,E31), (E12,E22,E23), (E31,E23,E33))))

# Transversely isotropic directional enhancements factors

m   = m3 # Symmetry axis

Emm = Constant(E33) # compressional enhancement factor
Emt = Constant(E31) # shear enhancement factor

#------------------
# Mesh, function spaces, BCs
#------------------

N = 5
mesh = UnitCubeMesh(N,N,N)

# Function spaces

U = VectorElement("Lagrange", mesh.ufl_cell(), 2) 
P = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, U*P) 

(u0, p0) = TrialFunctions(W) # the unknowns 
(v, q)   = TestFunctions(W) # the weight functions

w_lin  = Function(W)
w      = Function(W)
(u, p) = split(w) # the unknowns for NONLINEAR solver

# Boundary conditions (no slip at bottom face, free side faces, forced top face)

idx = 2 # vertical direction
def bnd_bot(x, on_boundary): return on_boundary and near(x[idx], 0) 
def bnd_top(x, on_boundary): return on_boundary and near(x[idx], 1) 

bc = [DirichletBC(W.sub(0), Constant((0,0,0)), bnd_bot)]
if dtype=='shear': bc += [DirichletBC(W.sub(0).sub(idx), Constant((0)), bnd_top)]

class Bnd_top(SubDomain): 
    def inside(self, x, on_boundary): return on_boundary and near(x[idx],1)

# Update "ds"

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0) 
boundaries.set_all(0)
Top = Bnd_top()
bnd_top_id = 1
Top.mark(boundaries, bnd_top_id)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#------------------
# Flow laws
#------------------

def eps(u):     
    return sym(grad(u))

def Aglen_factor(nglen): 
    return 2**((1-nglen)/(2*nglen)) # A --> Aglen_factor*A  if table values of A are used.

### Isotropic (Glen's law)

def viscosity_isotropic(u, A,n): 
    
    return Aglen_factor(n) * A**(-1/n) * inner(eps(u),eps(u))**((1-n)/(2*n)) 

def flowlaw_isotropic(u, A,n):   
    
    return viscosity_isotropic(u,A,n) * eps(u)

### Transversely isotropic (Johnson's law)

def viscosity_tranisotropic(u, A,n, M,Emt,Emm): 
    
    (kI,kM,kL) = tranisotropic_coefs(n, Emt, Emm)
    
    eps_ = eps(u)
    I2_  = inner(eps_,eps_)
    I4_  = inner(eps_,M)
    I5_  = inner(dot(eps_,eps_),M) 
    
    return Aglen_factor(n) * A**(-1/n) * (I2_ + kM * I4_**2 + 2*kL*I5_)**((1-n)/(2*n)) 
        
def flowlaw_tranisotropic(u, A,n, m,Emt,Emm): 
    
    (kI,kM,kL) = tranisotropic_coefs(n, Emt, Emm)
    
    eps_ = eps(u)
    M    = outer(m,m)
    I4_  = inner(eps_,M)
    
    FL  = eps_ + kI*I4_*Identity(3) + kM*I4_*M + kL*(dot(eps_,M)+dot(M,eps_)) 

    return viscosity_tranisotropic(u, A,n, M,Emt,Emm) * FL

def tranisotropic_coefs(n, Emt, Emm):
    
    Eexpo = 2.0/(n + 1.0) 
    
    kI = - (1/Emm**Eexpo-1)/2 
    kM = + (3*(1/Emm**Eexpo-1)-4*(1/Emt**Eexpo-1))/2 
    kL = + (1/Emt**Eexpo-1) 

    return (kI,kM,kL)

### Orthotropic (Gillet-Chaulet and others (2005), Rathmann and Lilien (2021))

def viscosity_orthotropic(u, A,n, J23,J31,J12,J4,J5,J6, Eij): 

    (lam1,lam2,lam3,lam4,lam5,lam6, gam) = orthotropic_coefs(n, Eij)
    
    return Aglen_factor(n) * A**(-1/n) * (\
                + lam1/gam * J23**2 \
                + lam2/gam * J31**2 \
                + lam3/gam * J12**2 \
                + 4 * (1/lam4) * J4**2 \
                + 4 * (1/lam5) * J5**2 \
                + 4 * (1/lam6) * J6**2 \
           )**((1-n)/(2*n)) 
        
def flowlaw_orthotropic(u, A,n, m1,m2,m3, Eij): 
    
    (lam1,lam2,lam3,lam4,lam5,lam6, gam) = orthotropic_coefs(n, Eij)
    
    eps_ = eps(u)
    
    M11 = outer(m1,m1)
    M22 = outer(m2,m2)
    M33 = outer(m3,m3)
    M23 = outer(m2,m3)
    M31 = outer(m3,m1)
    M12 = outer(m1,m2)
    
    eps11 = inner(eps_, M11)
    eps22 = inner(eps_, M22)
    eps33 = inner(eps_, M33)
    eps23 = inner(eps_, M23)
    eps31 = inner(eps_, M31)
    eps12 = inner(eps_, M12)
   
    J1 = (eps22 - eps33)/2.0 
    J2 = (eps33 - eps11)/2.0 
    J3 = (eps11 - eps22)/2.0 
    
    # Strictly speaking, these are sqrt(I4), sqrt(I5), and sqrt(I6)
    J4 = eps23
    J5 = eps31
    J6 = eps12
    
    J23 = -3/2.0 * eps11 # J23 := J2-J3 = (eps22 + eps33 - 2*eps11)/2.0d0 = -3/2.0d0 * eps11
    J31 = -3/2.0 * eps22 # J31 := J3-J1 = (eps11 + eps33 - 2*eps22)/2.0d0 = -3/2.0d0 * eps22
    J12 = -3/2.0 * eps33 # J12 := J1-J2 = (eps11 + eps22 - 2*eps33)/2.0d0 = -3/2.0d0 * eps33
    
    FL  = \
        + lam1/gam * J23 * (Identity(3) - 3*M11)/2 \
        + lam2/gam * J31 * (Identity(3) - 3*M22)/2 \
        + lam3/gam * J12 * (Identity(3) - 3*M33)/2 \
        + 4 * (1/lam4) * J4 * sym(M23) \
        + 4 * (1/lam5) * J5 * sym(M31) \
        + 4 * (1/lam6) * J6 * sym(M12)  

    return viscosity_orthotropic(u, A,n, J23,J31,J12,J4,J5,J6, Eij) * FL

def orthotropic_coefs(n, Eij):
    
    Eexpo = 2.0/(n + 1.0) 
    
    E11n = Eij[0,0]**Eexpo
    E22n = Eij[1,1]**Eexpo
    E33n = Eij[2,2]**Eexpo

    lam1 = 4/3.0*(-E11n+E22n+E33n)
    lam2 = 4/3.0*(+E11n-E22n+E33n)
    lam3 = 4/3.0*(+E11n+E22n-E33n)
    lam4 = 2* Eij[1,2]**Eexpo
    lam5 = 2* Eij[2,0]**Eexpo
    lam6 = 2* Eij[0,1]**Eexpo
    
    gam = -(E11n**2 + E22n**2 + E33n**2) + 2*(E11n*E22n + E11n*E33n + E22n*E33n)

    return (lam1,lam2,lam3,lam4,lam5,lam6, gam)

#------------------
# Weak form
#------------------
    
def weak_form(flowlaw, u,v, p,q, Aglen,nglen):
    
    a = (-div(v)*p + q*div(u) )*dx  
    L = rhoi*inner(g, v)*dx 
    
    if flowlaw == 'iso':   a += inner(    flowlaw_isotropic(u,Aglen,nglen),eps(v))*dx  
    if flowlaw == 'tiso':  a += inner(flowlaw_tranisotropic(u,Aglen,nglen, m,Emt,Emm),eps(v))*dx  
    if flowlaw == 'ortho': a += inner(  flowlaw_orthotropic(u,Aglen,nglen, m1,m2,m3, Eij),eps(v))*dx  

    L += dot(v,Constant(taushearmax)*r)*ds(bnd_top_id) # Traction on top face 
            
    return (a,L)

def solve_Stokes(flowlaw, w_guess=None, tolfac=1, relaxation=0.5, maxiter=50):
   
    # Weak Stokes problem
    (a_lin,L_lin) = weak_form(flowlaw, u0,v, p0,q, Aglen, 1) # nglen = 1
    solve(a_lin==L_lin, w_lin, bc) # Linear solve
    w.assign(w_lin)

    # If nonlinear solve, then use the init guess stored in "w"
    if nglen > 1:    
        (a_nlin,L_nlin) = weak_form(flowlaw, u,v, p,q, Aglen,nglen)
        F = a_nlin-L_nlin
        solve(F == 0, w, bc, solver_parameters={"newton_solver":  {"relative_tolerance": tolfac*1e-4,"absolute_tolerance": tolfac*1e-3,"relaxation_parameter":relaxation, "maximum_iterations":maxiter}})
        
    return w

#------------------
# Cube deformation experiments using different rheologies
#------------------

def plot_2d_slice(mesh, u_iso, u_tiso, u_ortho):
    
    bmesh = BoundaryMesh(mesh, "exterior") # Create a UnitSquareMesh of topology 2 in 3 dimensions
    ff = MeshFunction("size_t", bmesh, 2)
    ff.set_all(0)
    y0 = AutoSubDomain(lambda x, on_bnd: near(x[1],0)) # Choose a side with normal in y-direction
    y0.mark(ff, 1)
    smesh = SubMesh(bmesh, ff, 1)
    xy = smesh.coordinates()
    xy[:, 1] = y0_plane # Move the 2D slice
        
    Q = FunctionSpace(smesh, 'CG', 2) # interpolate from u0 to the slice
    
    ux_iso, uy_iso, uz_iso       = Function(Q), Function(Q), Function(Q)
    ux_tiso, uy_tiso, uz_tiso    = Function(Q), Function(Q), Function(Q)
    ux_ortho, uy_ortho, uz_ortho = Function(Q), Function(Q), Function(Q)
    
    LagrangeInterpolator.interpolate(ux_iso, u_iso.sub(0))
    LagrangeInterpolator.interpolate(uy_iso, u_iso.sub(1))
    LagrangeInterpolator.interpolate(uz_iso, u_iso.sub(2))
    
    LagrangeInterpolator.interpolate(ux_tiso, u_tiso.sub(0))
    LagrangeInterpolator.interpolate(uy_tiso, u_tiso.sub(1))
    LagrangeInterpolator.interpolate(uz_tiso, u_tiso.sub(2))
    
    LagrangeInterpolator.interpolate(ux_ortho, u_ortho.sub(0))
    LagrangeInterpolator.interpolate(uy_ortho, u_ortho.sub(1))
    LagrangeInterpolator.interpolate(uz_ortho, u_ortho.sub(2))
    
    qq = (dtype=='comp')
    plt = vplot(ux_iso, at=0, shape=(3,3), cmap='RdBu' if qq else 'Blues', elevation=-90, interactive = False)
    plt = vplot(uy_iso, at=1, shape=(3,3), cmap='RdBu', interactive = False)
    plt = vplot(uz_iso, at=2, shape=(3,3), cmap='Reds_r' if qq else 'RdBu', interactive = False)

    plt = vplot(ux_tiso, at=3, shape=(3,3), cmap='RdBu' if qq else 'Blues', interactive = False)
    plt = vplot(uy_tiso, at=4, shape=(3,3), cmap='RdBu', interactive = False)
    plt = vplot(uz_tiso, at=5, shape=(3,3), cmap='Reds_r' if qq else 'RdBu', interactive = False)
    
    plt = vplot(ux_ortho, at=6, shape=(3,3), cmap='RdBu' if qq else 'Blues', interactive = False)
    plt = vplot(uy_ortho, at=7, shape=(3,3), cmap='RdBu', interactive = False)
    plt = vplot(uz_ortho, at=8, shape=(3,3), cmap='Reds_r' if qq else 'RdBu', interactive = False)
    
    plt.screenshot('%s_nglen%i.png'%(dtype, nglen))

# ---------- MAIN ------------

w_iso = solve_Stokes('iso')
(u_iso,_) = w_iso.split()

w_tiso = solve_Stokes('tiso')
(u_tiso,_) = w_tiso.split()

w_ortho = solve_Stokes('ortho')
(u_ortho,_) = w_ortho.split()

plot_2d_slice(mesh, u_iso, u_tiso, u_ortho)
