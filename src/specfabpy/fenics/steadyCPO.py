#!/usr/bin/python3
# Nicholas Rathmann <rathmann@nbi.ku.dk>, 2024-

"""
Steady SSA CPO solver for ice
"""

import os, sys, copy, code # code.interact(local=locals())
import numpy as np
from tabulate import tabulate
import copy, xarray, pickle, pyproj # pyproj needed!
from dolfin import *
from .ice import IceFabric
from ..specfabpy import specfabpy as sf__ # sf private copy 
import matplotlib.tri as tri

class SteadyCPO():

    uxname, uyname = 'VX', 'VY'
    xuname, yuname = 'x', 'y'

    CAFFE_params = (0.1, 10) # Emin, Emax of CAFFE

    velscale = 3.17098e+8 # m/s to m/yr
    lenscale = 1e-3 # m to km

    modelplane = 'xz' # xz is much faster than xy (renders state vector real-valued)

    bc_isotropic  = 0
    bc_zsinglemax = 1

    def __init__(self, domain):
    
        self.exname = domain['name'] # experiment name
        self.fvel   = domain['fmeasures']
        self.fbed   = domain['fbedmachine']
        self.fmesh  = domain['fgeo'][:-4]
        self.dxy_u  = domain['subsample_u']
        self.dxy_h  = domain['subsample_h']
        
        ### Determine bounding box of domain for interpolation

        content = open(domain['fgeo']).readlines()
        p0 = (self._getnum(content[2]), self._getnum(content[3]))
        p1 = (self._getnum(content[4]), self._getnum(content[5]))
        #print(p0,p1)
        # ...unpack
        self.x0, self.y0 = p0
        self.x1, self.y1 = p1
        self.dx = 0.025 * abs(self.x1-self.x0)
        self.dy = 0.025 * abs(self.y1-self.y0)
        
    def _getnum(self, s):
        part1, part2 = s.split("=", 1) # Split the string into two parts based on the start delimiter
        result = part2.split(";", 1)[0] # Split the second part based on the end delimiter
        return float(result)
        
    def preprocess(self):

        """
        Slice and subsample maps of velocity and geometry.
        Results are saved to .h5 and .pkl files.
        """
    
        print('*** Preprocessing started')
        print(tabulate([['name', self.exname], ['fmeasures', self.fvel], ['fbedmachine', self.fbed], \
                        ['subsample_u', self.dxy_u], ['subsample_h', self.dxy_h]], \
                        headers=['Parameter', 'Value'], tablefmt="rst"))

        ### Make mesh

        print('*** Generating finite element mesh')
        os.system('gmsh -2 -format msh2 %s.geo'%(self.fmesh))
        os.system('dolfin-convert %s.msh %s.xml'%(self.fmesh, self.fmesh))

        ### Load mesh
        
        mesh,boundaries, Q,Q2,V, coords,cells = self.get_mesh()

        # xarray objects for interpolation points (used below)
        coordsQ = Q.tabulate_dof_coordinates().reshape((-1, 2)).T
        xQ = xarray.DataArray(coordsQ[0,:], dims="z")
        yQ = xarray.DataArray(coordsQ[1,:], dims="z")
        
        coordsQ2 = Q2.tabulate_dof_coordinates().reshape((-1, 2)).T
        xQ2 = xarray.DataArray(coordsQ2[0,:], dims="z")
        yQ2 = xarray.DataArray(coordsQ2[1,:], dims="z")

        ### Slice and coarsen

        print('*** Coarsening input fields')

        ds_vel0 = xarray.open_mfdataset(self.fvel)
        ds_bed0 = xarray.open_mfdataset(self.fbed)
        xrng = slice(self.x0-self.dx, self.x1+self.dx)
        yrng = slice(self.y1+self.dy, self.y0-self.dy)
        ds_u = ds_vel0.sel(x=xrng,y=yrng).coarsen(x=self.dxy_u, boundary='trim').mean().coarsen(y=self.dxy_u, boundary='trim').mean()
        ds_b = ds_bed0.sel(x=xrng,y=yrng).coarsen(x=self.dxy_h, boundary='trim').mean().coarsen(y=self.dxy_h, boundary='trim').mean()

        ### Geometry
        
        print('*** Processing geometry')
        
        S, B, H, mask = Function(Q), Function(Q), Function(Q), Function(Q)
        kwargsLN = dict(x=xQ, y=yQ, method='linear')
        kwargsNN = dict(x=xQ, y=yQ, method='nearest')
        S.vector()[:]    = ds_b.surface.interp(**kwargsLN).to_numpy()
        B.vector()[:]    = ds_b.bed.interp(**kwargsLN).to_numpy()
        H.vector()[:]    = ds_b.thickness.interp(**kwargsLN).to_numpy()
        mask.vector()[:] = ds_b.mask.interp(**kwargsNN).to_numpy()

        ### Velocity
        
        print('*** Processing velocity')

        ds_ui = ds_u.interp(x=xQ2, y=yQ2, method='linear')
        ux_np = 1/self.velscale * getattr(ds_ui, self.uxname).to_numpy()
        uy_np = 1/self.velscale * getattr(ds_ui, self.uyname).to_numpy()
        I = np.argwhere(np.isnan(ux_np))
        if len(I) > 0:
            print('*** Warning: NaNs found in velocity map where nodes are located... please make sure the model domain does not extend outside that of the velocity product')
            if 1:
                print('...setting %i components to zero and proceeding silently.'%(len(I)))
                ux_np[I] = 1e-15
                uy_np[I] = 1e-15
            else:
                print('... exiting. I=', I.T, coordsQ2[:,I])
                sys.exit(1)
        ux, uy, u = Function(Q2), Function(Q2), Function(V)
        ux.vector()[:] = ux_np
        uy.vector()[:] = uy_np
        FunctionAssigner(V, [Q2,Q2]).assign(u, [ux,uy]) # set vector field components

        umag = Function(Q2)
        umag_np = np.sqrt(np.power(ux_np,2)+np.power(uy_np,2))
        umag.vector()[:] = umag_np
    
        D    = project(sym(grad(u)), TensorFunctionSpace(mesh, 'CG', 1))
        epsE = project(inner(D,D)+tr(D)**2, Q2)
        epsE.vector()[:] = np.sqrt(epsE.vector()[:]/2) # sqrt(D:D/2)

        if np.isnan(umag.vector()[:]).any(): raise ValueError('*** Warning: umag contains NaNs')
        if (umag.vector()[:]<0).any():       raise ValueError('*** Warning: umag contains values <0 ')

        if np.isnan(epsE.vector()[:]).any(): raise ValueError('*** Warning: epsE contains NaNs')
        if (epsE.vector()[:]<0).any():       raise ValueError('*** Warning: epsE contains values <0 ')

        ### Save FEM fields

        print('*** Saving FEM fields to %s-inputs.h5'%(self.exname))

        fout = '%s-inputs.h5'%(self.exname)
        f = HDF5File(MPI.comm_world, fout, "w")
        f.write(u, "/u")
        f.write(umag, "/umag")
        f.write(epsE, "/epsE")
        f.write(S, "/S")
        f.write(B, "/B")
        f.write(H, "/H")
        f.write(mask, "/mask")
        f.close()
        
        ### Save numpy fields

        print('*** Saving numpy fields to %s-inputs.pkl'%(self.exname))

        fout = '%s-inputs.pkl'%(self.exname)
        with open(fout, 'wb') as handle:
            data = [coords,cells] + [F.compute_vertex_values(mesh) for F in (ux,uy,umag,epsE, S,B,H,mask)]
            pickle.dump(data, handle)

    def get_mesh(self, mapscale=1):
    
        mesh = Mesh('%s.xml'%(self.fmesh)) 
        boundaries = MeshFunction('size_t', mesh, '%s_facet_region.xml'%(self.fmesh))
        Q  = FunctionSpace(mesh, 'CG', 1)
        Q2 = FunctionSpace(mesh, 'CG', 2) # this space is used when constructing the FEM velocity field (strain-rate tensor requires it to be DG2; Gibbs phenomenon?) 
        V  = VectorFunctionSpace(mesh, 'CG', 2)
        coords = mapscale*copy.deepcopy(mesh.coordinates().reshape((-1, 2)).T)
        cells = mesh.cells()
        return (mesh, boundaries, Q,Q2,V, coords,cells)
        
    def triang(self, coords, cells, mapscale=1):
    
        return tri.Triangulation(mapscale*coords[0], mapscale*coords[1], triangles=cells)
        
    def feminputs(self, Q=None, Q2=None, V=None):
    
        if Q is None: mesh,boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u,umag,epsE = Function(V), Function(Q2), Function(Q2)
        S,B,H,mask = Function(Q), Function(Q), Function(Q), Function(Q)
        f = HDF5File(MPI.comm_world,'%s-inputs.h5'%(self.exname),'r')
        f.read(u, '/u')
        f.read(umag, '/umag')
        f.read(epsE, '/epsE')
        f.read(S, '/S')
        f.read(B, '/B')
        f.read(H, '/H')
        f.read(mask, '/mask')
        f.close()
        return (u,umag,epsE, S,B,H,mask)
        
    def npinputs(self):
    
        #print('*** Loading numpy inputs')
        with open('%s-inputs.pkl'%(self.exname), 'rb') as handle:
            return pickle.load(handle)
        
    """
    Solver
    """

    def solve(self, problem, numerics):
    
        isDDRX = len(np.shape(problem['T']))>0
        print(tabulate([['PROBLEM', 'LROT+ADVEC' if not isDDRX else 'LROT+DDRX+ADVEC'], 
                        ['T',          'None' if not isDDRX else ','.join(['%.1f'%(_) for _ in problem['T']])], 
                        ['L',          numerics['L']], 
                        ['nu_orimul',  '%.2e'%(numerics['nu_orimul'])], 
                        ['nu_real',    '%.2e'%(numerics['nu_real'])], 
                        ['nu_realmul', '1' if not isDDRX else ','.join(['%.1f'%(_) for _ in numerics['nu_realmul']])], 
                        ['modelplane', self.modelplane], 
                        ], headers=['Parameter', 'Value'], tablefmt="rst"))
    
        ### Mesh

        print('*** Loading mesh and input fields')
    
        mesh, boundaries, Q,Q2,V, coords,cells = self.get_mesh()
        u, umag, *_ = self.feminputs(Q=Q, Q2=Q2, V=V)        
        
        ### Init ice fabric class
        
        fab = IceFabric(mesh, boundaries, numerics['L'], nu_realspace=numerics['nu_real'], nu_multiplier=numerics['nu_orimul'], \
                            modelplane=self.modelplane, CAFFE_params=self.CAFFE_params)
        
        ### Boundary conditions
        
        nlm_list = [np.array(fab.nlm_iso), fab.sf.nlm_ideal([0,1,0], 0, numerics['L'])] # isotropic or z-SMAX (note y is vertical since modelplane is xz)
        bc_vals  = [np.real(fab.sf.nlm_to_rnlm(nlm_list[bc[1]], fab.sf.get_rnlm_len())) for bc in problem['bcs']]
        bc_ids   = [bc[0] for bc in problem['bcs']]

        fab.set_BCs(bc_vals, [0*val for val in bc_vals], bc_ids, domain=boundaries) # (real, imag) parts on bc_ids
                
        ### Solve steady SSA problem
        
        #print('*** Solving for steady SSA CPO field')
        
        if not isDDRX: Gamma0 = None
        else:          Gamma0 = [fab.Gamma0_Lilien23_lab(u, _+273.15) for _ in problem['T']] # list of DDRX rate factors used to gradually approach solution 

        S = sym(grad(u)) # approximation strain rate tensor (D) as coaxial to stress tensor (S)
        fab.solvesteady(u, S, iota=+1, Gamma0=Gamma0, nu_realspacemul=numerics['nu_realmul'], LROT_guess=True)
        
        ### Calculate fabric-derived quantities
        
        mi, lami = fab.eigenframe(*coords) # sorted such that last entry is the out-out-model-plane eigenpair
        E = fab.get_E_CAFFE(u).compute_vertex_values(mesh)    
        E[np.isnan(E)] = 1 # not sure why a few values sometimes are nan
            
        ### Save numpy solution

        fout = '%s-solution-%s.pkl'%(self.exname, problem['name'])
        print('*** Saving numpy solution to %s'%(fout))
        with open(fout, 'wb') as handle:
            data = (coords,cells, mi,lami,E)
            pickle.dump(data, handle)
         
        ### Save FEM solution
        
        fout = '%s-solution-%s.h5'%(self.exname, problem['name'])
        print('*** Saving FEM solution to %s'%(fout))
        f = HDF5File(MPI.comm_world, fout, "w")
        f.write(fab.s,"/s")
        f.close()
        
    def npsolution(self, probname):
    
        fout = '%s-solution-%s.pkl'%(self.exname, probname)
        with open(fout, 'rb') as handle:
            return pickle.load(handle)
       
    def femsolution(self, probname, mesh=None, boundaries=None):
        ### NOT USED, DELETE?
        
        if mesh is None: mesh, boundaries, *_ = self.get_mesh()
        fab = IceFabric(mesh, boundaries, self.L)
        s = Function(fab.S)
        fout = '%s-solution-%s.h5'%(self.exname, probname)
        f = HDF5File(MPI.comm_world,fout,"r")
        f.read(s,"/s")
        f.close()
        fab.s.assign(s)
        return fab
            
    """
    AUX
    """

    def bmesh(self, bcs, mapscale=1):
    
        mesh, boundaries, *_ = self.get_mesh()
        boundarymesh = BoundaryMesh(mesh, 'exterior')
        bdim = boundarymesh.topology().dim()
        boundary_boundaries = MeshFunction('size_t', boundarymesh, bdim)
        boundary_boundaries.set_all(0)
        for i, facet in enumerate(entities(boundarymesh, bdim)):
            parent_meshentity = boundarymesh.entity_map(bdim)[i]
            parent_boundarynumber = boundaries.array()[parent_meshentity]
            boundary_boundaries.array()[i] = parent_boundarynumber
        bmeshes = [SubMesh(boundarymesh, boundary_boundaries, bc[0]) for bc in bcs]
        coords = [copy.deepcopy(bmesh.coordinates().reshape((-1, 2)).T) for bmesh in bmeshes]
        return (coords, bmeshes)
        
    def xyboundaries(self):
    
        ### Not guaranteed to give the correct connectivity between nodes on the boundary (use with caution)
        (coords, bmeshes) = self.bmesh()
        x, y = [], []
        for c in coords:
            xy = c.T
            tour, total_distance = self.tsp(xy)
            xi = np.array([xy[i][0]  for i in tour])
            yi = np.array([xy[i][1]  for i in tour])
            x.append(xi)
            y.append(yi)
        return (x,y)
        
    """
    TSP for determining ordered list of boundary coordinates 
    Code from https://www.askpython.com/python/examples/travelling-salesman-problem-python
    """
    
    def distance(self, city1, city2):
      # Replace this with your distance calculation function (e.g., Euclidean distance)
      x1, y1 = city1
      x2, y2 = city2
      return ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
     
    def tsp(self, cities):
      visited = [False] * len(cities)
      current_city = 0
     
      tour = []
      total_distance = 0
     
      for _ in range(len(cities)):
        visited[current_city] = True
        tour.append(current_city)
     
        next_city = None
        min_distance = float('inf')  # Represents positive infinity
     
        for i in range(len(cities)):
          if visited[i]:
            continue
     
          d = self.distance(cities[current_city], cities[i])
          if d < min_distance:
            min_distance = d
            next_city = i
     
        current_city = next_city
        total_distance += min_distance
     
      return tour, total_distance

