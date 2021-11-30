# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2019-2021

#-------------------
# Core
#-------------------

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC
OPTS=-O2 -mcmodel=medium -lm -llapack -lblas -lfftw3
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts
RHEOLOGIES=rheologies

PATHDEMO=../demo
PATHFABEVO=$(PATHDEMO)/fabric-evolution
PATHRHEOTEST=$(PATHDEMO)/rheology-test

ALLOBJS=$(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(GAUNT).o $(RHEOLOGIES).o
ALLMODS=$(ALLOBJS:.o=.mod)
ALLSRCS=$(ALLOBJS:.o=.f90)

#-------------------
# Elmer, JOSEF, ...
#-------------------

ELMERSRCS=$(wildcard $(PATHDEMO)/elmer/elmer*.f90)
ELMERS=$(ELMERSRCS:.f90=)

########################################

#-------------------
# Core
#-------------------

fabric-evolution-demo: $(SPECFAB).o
	$(COMPILER) $(PATHFABEVO)/run-fabric-evolution-latrot.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(PATHFABEVO)/run-fabric-evolution-latrot 
	$(COMPILER) $(PATHFABEVO)/run-fabric-evolution-ddrx.f90   $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(PATHFABEVO)/run-fabric-evolution-ddrx
	mkdir -p $(PATHFABEVO)/solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-latrot  uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-latrot  ss_xz :: simple shear (ss) along the x-z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-ddrx uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-ddrx ss_xz :: simple shear (ss) along the x-z plane"

specfabpy: $(SPECFAB).o
	f2py -lm -llapack -lblas -lfftw3 -I. $(ALLOBJS) -c -m specfabpy specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet 
	cp specfabpy.cpython* $(PATHDEMO)/
	@echo "\n*** Python module installed in demo path."
	mkdir -p $(PATHFABEVO)/solutions
	@echo "\n*** To get going, try the demo:"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py ss_xz ::: for simple shear (ss) along the x-z plane"

lib$(SPECFAB).so: $(ALLSRCS) $(ALLOBJS)
	$(COMPILER) $(OPTS) -shared $(ALLSRCS) -o $@

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(RHEOLOGIES).f90
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

$(MOMENTS).o: 
	@echo "*** Compiling structure tensor expressions and gaunt coefficients... this may take some time but is required only once ***"
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(GAUNT).o: 
	$(COMPILER) $(OPTS) -c $(GAUNT).f90

clean:
	rm -f $(PATHFABEVO)/run-fabric-evolution-latrot $(PATHFABEVO)/run-fabric-evolution-ddrx *.o *.mod *.so
	
clear:
	rm -f $(PATHFABEVO)/run-fabric-evolution-latrot $(PATHFABEVO)/run-fabric-evolution-ddrx $(SPECFAB).o $(SPECFAB).mod *.so

rheology-test: $(SPECFAB).o
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-orthotropic.f90           $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-orthotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-transverselyisotropic.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-transverselyisotropic
	
#-------------------
# Elmer, JOSEF, ...
#-------------------

#elmer: $(ELMERS)

#$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(ALLMODS)
#	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@


