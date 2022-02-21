# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2019-2022

#-------------------
# Compiler
#-------------------

COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC -Wno-integer-division -Wno-unused-dummy-argument
OPTS=-O2 -mcmodel=medium -lm -llapack -lblas
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

#-------------------
# Core modules
#-------------------

PATHDEMO=../demo
PATHTEST=../tests
PATHFABEVO=$(PATHDEMO)/fabric-evolution
PATHRHEOTEST=$(PATHTEST)/rheology-test

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts
HOMOGEN=homogenizations
ENHFAC=enhancementfactors
MANDEL=mandel
DYNAMICS=dynamics
RHEOLOGIES=rheologies
REDUCED=reducedform
LAMBDASOLVER=lambdasolver

OBJCORE=$(SPECFAB).o $(TENPROD).o $(HOMOGEN).o $(ENHFAC).o $(MANDEL).o $(MOMENTS).o $(GAUNT).o $(DYNAMICS).o $(RHEOLOGIES).o $(REDUCED).o $(LAMBDASOLVER).o
MODCORE=$(OBJCORE:.o=.mod)
SRCCORE=$(OBJCORE:.o=.f90)

#-------------------
# Math77 library
#-------------------

PATHM77=include/math77
OBJM77=$(PATHM77)/amach.o $(PATHM77)/derv1.o $(PATHM77)/dnqsol.o $(PATHM77)/dnrm2.o $(PATHM77)/erfin.o $(PATHM77)/ermsg.o $(PATHM77)/ierm1.o $(PATHM77)/ierv1.o
SRCM77=$(OBJM77:.o=.f)

#-------------------
# Elmer, JOSEF, ...
#-------------------

ELMERSRCS=$(wildcard $(PATHDEMO)/elmer/elmer*.f90)
ELMERS=$(ELMERSRCS:.f90=)

#-------------------
# Collect all
#-------------------

OBJALL=$(OBJCORE) $(OBJM77)
MODALL=$(MODCORE)
SRCALL=$(SRCCORE) $(SRCM77)

########################################

#-------------------
# Make rules
#-------------------

fabric-evolution-demo: $(SPECFAB).o
	$(COMPILER) $(PATHFABEVO)/run-fabric-evolution-latrot.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHFABEVO)/run-fabric-evolution-latrot 
	$(COMPILER) $(PATHFABEVO)/run-fabric-evolution-ddrx.f90   $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHFABEVO)/run-fabric-evolution-ddrx
	@mkdir -p $(PATHFABEVO)/solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-latrot  uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-latrot  ss_xz :: simple shear (ss) along the x-z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-ddrx uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); ./run-fabric-evolution-ddrx ss_xz :: simple shear (ss) along the x-z plane"

specfabpy: $(SPECFAB).o
	rm -f specfabpy.pyf
	f2py --no-lower -m specfabpy -h specfabpy.pyf specfabpy.f90 --quiet 
	f2py -lm -llapack -lblas -I. $(OBJALL) -c specfabpy.pyf specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet
	@cp specfabpy.cpython* $(PATHDEMO)/
	@echo "\n*** Python module (specfabpy.cpython) installed in demo path."
	@mkdir -p $(PATHFABEVO)/solutions
	@echo "\n*** To get going, try the demo:"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py ss_xz ::: for simple shear (ss) along the x-z plane"

rheology-test: $(SPECFAB).o
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-orthotropic.f90           $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-orthotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-transverselyisotropic.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-transverselyisotropic

aux-tests: $(SPECFAB).o
	$(COMPILER) $(PATHTEST)/reduced-coef-test/reduced-coef-test.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHTEST)/reduced-coef-test/reduced-coef-test

lib$(SPECFAB).so: $(SRCALL) $(OBJALL)
	$(COMPILER) $(OPTS) -shared $(SRCALL) -o $@

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o $(LAMBDASOLVER).o
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(RHEOLOGIES).f90 
	$(COMPILER) $(OPTS) -c $(HOMOGEN).f90
	$(COMPILER) $(OPTS) -c $(ENHFAC).f90
	$(COMPILER) $(OPTS) -c $(DYNAMICS).f90
	$(COMPILER) $(OPTS) -c $(REDUCED).f90
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

$(MOMENTS).o: 
	@echo "*** Compiling structure tensor expressions and gaunt coefficients... this may take some time but is required only once ***"
	$(COMPILER) $(OPTS) -c $(MANDEL).f90 
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(GAUNT).o: 
	$(COMPILER) $(OPTS) -c $(GAUNT).f90

$(LAMBDASOLVER).o:
	rm -f $(OBJM77) $(LAMBDASOLVER).o $(LAMBDASOLVER).mod $(PATHM77)/test
	cd $(PATHM77) && $(COMPILER) $(OPTS) -c *.f
	$(COMPILER) $(OPTS) -c $(LAMBDASOLVER).f90
	$(COMPILER) $(OPTS) $(PATHM77)/test.f90 $(LAMBDASOLVER).o $(OBJM77) -o $(PATHM77)/test

clear:
	rm -f $(PATHRHEOTEST)/rheology-test-orthotropic $(PATHRHEOTEST)/rheology-test-transverselyisotropic
	rm -f $(PATHFABEVO)/run-fabric-evolution-latrot $(PATHFABEVO)/run-fabric-evolution-ddrx 
	rm -f $(SPECFAB).o $(SPECFAB).mod *.so
	rm -f $(LAMBDASOLVER).o $(LAMBDASOLVER).mod $(OBJM77) $(PATHM77)/test

clean: clear
	rm -f *.o *.mod *.so

#-------------------
# Elmer, JOSEF, ...
#-------------------

#elmer: $(ELMERS)

#$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(MODCORE)
#	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@


