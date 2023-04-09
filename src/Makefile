# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

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
PATHRHEOTEST=$(PATHTEST)/constitutive-equations

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts
HOMOGEN=homogenizations
ENHFAC=enhancementfactors
MANDEL=mandel
DYNAMICS=dynamics
RHEOLOGIES=rheologies
ELASTICITIES=elasticities
WAVEPROP=wavepropagation
ROTATION=rotation
REDUCED=reducedform
LAMBDASOLVER=lambdasolver

DEPLESS=$(TENPROD).o $(MANDEL).o $(REDUCED).o $(GAUNT).o
WITHDEPS=$(RHEOLOGIES).o $(ELASTICITIES).o $(MOMENTS).o $(LAMBDASOLVER).o $(DYNAMICS).o $(HOMOGEN).o $(ENHFAC).o $(ROTATION).o $(WAVEPROP).o
# 
OBJCORE=$(DEPLESS) $(WITHDEPS)
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

ELMERSRCS=$(wildcard elmer/specfab*.f90)
ELMERS=$(ELMERSRCS:.f90=)

#-------------------
# Collect all
#-------------------

OBJALL=$(SPECFAB).o $(OBJCORE) $(OBJM77)
SRCALL=$(SPECFAB).f90 $(SRCCORE) $(SRCM77)

#-------------------
# specfabpy Python versions
#-------------------

#PYVERS=3.6 3.7 3.8 3.9 3.10 3.11
PYVERS=3

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
	for pyver in $(PYVERS); do \
		echo "*** Making specfabpy for python$$pyver"; \
		rm -f specfabpy.pyf; \
		python$$pyver -m numpy.f2py --no-lower -m specfabpy -h specfabpy.pyf specfabpy.f90 --quiet; \
		python$$pyver -m numpy.f2py -lm -llapack -lblas -I. $(OBJALL) -c specfabpy.pyf specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet; \
	done
	@cp specfabpy.cpython* $(PATHDEMO)/
	@echo "\n*** Python module (specfabpy) installed in /demo"
	@mkdir -p $(PATHFABEVO)/solutions
	@echo "\n*** To get going, try the demo:"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(PATHFABEVO); python3 run-fabric-evolution-latrot.py ss_xz ::: for simple shear (ss) along the x-z plane"

python: specfabpy
	python setup.py install

constitutive-tests: $(SPECFAB).o
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-orthotropic.f90             $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-orthotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-transverselyisotropic.f90   $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-transverselyisotropic
	$(COMPILER) $(PATHRHEOTEST)/elasticity-test-transverselyisotropic.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/elasticity-test-transverselyisotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-isotropic-porous.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-isotropic-porous

aux-tests: $(SPECFAB).o
	$(COMPILER) $(PATHTEST)/reduced-coef-test/reduced-coef-test.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHTEST)/reduced-coef-test/reduced-coef-test
	$(COMPILER) $(PATHTEST)/wigner-d-rotation-test/discrete-rotation-test.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHTEST)/wigner-d-rotation-test/discrete-rotation-test

lib$(SPECFAB).so: $(SRCALL) $(OBJALL)
	$(COMPILER) $(OPTS) -shared $(SRCALL) -o $@

$(SPECFAB).o: $(SPECFAB).f90 $(OBJCORE) elmer/specfab_elmer.f90 josef/specfab_josef.f90
	$(COMPILER) $(OPTS) -c $<

$(DEPLESS): %.o: %.f90
	$(COMPILER) $(OPTS) -c $<

$(RHEOLOGIES).o: $(RHEOLOGIES).f90 $(TENPROD).o $(LAMBDASOLVER).o
	$(COMPILER) $(OPTS) -c $<

$(ELASTICITIES).o: $(ELASTICITIES).f90 $(TENPROD).o
	$(COMPILER) $(OPTS) -c $<

$(DYNAMICS).o: $(DYNAMICS).f90 tensorialdynamics.f90 $(TENPROD).o $(MOMENTS).o $(GAUNT).o
	$(COMPILER) $(OPTS) -c $<

$(HOMOGEN).o: $(HOMOGEN).f90 $(MOMENTS).o $(RHEOLOGIES).o $(ELASTICITIES).o
	$(COMPILER) $(OPTS) -c $<

$(ENHFAC).o: $(ENHFAC).f90 $(MOMENTS).o $(TENPROD).o $(HOMOGEN).o
	$(COMPILER) $(OPTS) -c $<

$(MOMENTS).o: $(MOMENTS).f90 $(MANDEL).o
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(WAVEPROP).o: $(WAVEPROP).f90 $(MOMENTS).o $(MANDEL).o $(TENPROD).o $(ELASTICITIES).o $(HOMOGEN).o
	$(COMPILER) $(OPTS) -c $<

$(ROTATION).o: $(ROTATION).f90 $(MOMENTS).o 
	$(COMPILER) $(OPTS) -c $<

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

all: $(SPECFAB).o lib$(SPECFAB).so

#-------------------
# Elmer, JOSEF, ...
#-------------------

#elmer: $(ELMERS)

#$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(MODCORE)
#	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@



