# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-

### Compiler

COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC -Wno-integer-division -Wno-unused-dummy-argument
# -mcmodel=medium might be needed if support for L>20 truncation is requested
OPTS=-O2 -mcmodel=small -lm -llapack -lblas
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

F2PY_VERBOSE=--quiet

### Core modules

PATHSPECFABPY=specfabpy
PATHDEMO=../demo
PATHTEST=../tests
PATHFABEVO=$(PATHDEMO)/fabric-evolution
PATHRHEOTEST=$(PATHTEST)/constitutive-equations

SPECFAB=specfab
HEADER=header
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
DEFORMMODES=deformationmodes
IDEALSTATE=idealstate
FRAMES=frames
LAMBDASOLVER=lambdasolver
GOLF=golf
ELMER=specfab_elmer

DEPLESS=$(HEADER).o $(TENPROD).o $(MANDEL).o $(REDUCED).o $(GAUNT).o $(GOLF).o
WITHDEPS=$(RHEOLOGIES).o $(ELASTICITIES).o $(MOMENTS).o $(LAMBDASOLVER).o $(DYNAMICS).o $(HOMOGEN).o $(ENHFAC).o $(ROTATION).o $(WAVEPROP).o $(DEFORMMODES).o $(IDEALSTATE).o $(FRAMES).o $(ELMER).o 

OBJCORE=$(DEPLESS) $(WITHDEPS)
SRCCORE=$(OBJCORE:.o=.f90)

### Math77 library

PATHM77=include/math77
OBJM77=$(PATHM77)/amach.o $(PATHM77)/derv1.o $(PATHM77)/dnqsol.o $(PATHM77)/dnrm2.o $(PATHM77)/erfin.o $(PATHM77)/ermsg.o $(PATHM77)/ierm1.o $(PATHM77)/ierv1.o
SRCM77=$(OBJM77:.o=.f)

### Elmer, ...

ELMERSRCS=$(wildcard specfab*_elmer.f90)
ELMERS=$(ELMERSRCS:.f90=)

### Collect all

OBJALL=$(SPECFAB).o $(OBJCORE) $(OBJM77)
SRCALL=$(SPECFAB).f90 $(SRCCORE) $(SRCM77)

### Python versions supported

PYVERS=3
#PYVERS=3.8 3.9 3.10 3.11 3.12 3.13

#---------------------------

all: specfabpy

# Linux and Mac
install: specfabpy
	python3 setup.py install --user
	# alternative on mac, run: pip3 install .

devinstall: specfabpy
	python3 setup-dev.py install --user

# Firedrake/icepack docker enviroment
icepack: specfabpy
	@echo "*** You can also just run: pip3 install specfabpy ***"
	#pip3 install setuptools==61 cartopy cmasher
	#pip3 install . 
	python3 setup.py install
	
specfabpy: $(SPECFAB).o
	@echo "*** Archiving objects into $(PWD)/libspecfab.a";
	rm -f $(PWD)/libspecfab.a;
	ar rcs $(PWD)/libspecfab.a $(OBJALL);
	for pyver in $(PYVERS); do \
		echo "*** Making specfabpy for python$$pyver"; \
		rm -f specfabpy.pyf; \
		python$$pyver -m numpy.f2py --no-lower -m specfabpy -h specfabpy.pyf specfabpy.f90 $(F2PY_VERBOSE); \
		python$$pyver -m numpy.f2py -lm -llapack -lblas -lspecfab -I$(PWD) -I$(PWD)/include -L$(PWD) $(OBJALL) -c specfabpy.pyf specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=small" $(F2PY_VERBOSE); \
	done
	@echo "*** Moving python *.so modules to src/$(PATHSPECFABPY)"
	@mv -f specfabpy.cpython* $(PATHSPECFABPY)/
	
clear:
	rm -f $(PATHRHEOTEST)/rheology-test-orthotropic $(PATHRHEOTEST)/rheology-test-transverselyisotropic
	rm -f $(PATHFABEVO)/run-fabric-evolution-latrot $(PATHFABEVO)/run-fabric-evolution-ddrx 
	rm -f $(SPECFAB).o $(SPECFAB).mod *.so $(PATHSPECFABPY)/*.so
	rm -f $(LAMBDASOLVER).o $(LAMBDASOLVER).mod $(OBJM77) $(PATHM77)/test $(FRAMES).o

clean: clear
	rm -f *.o *.mod *.so

# old, not used
#all: $(SPECFAB).o lib$(SPECFAB).so
	
#---------------------------

constitutive-tests: $(SPECFAB).o
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-orthotropic.f90             $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-orthotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-transverselyisotropic.f90   $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-transverselyisotropic
	$(COMPILER) $(PATHRHEOTEST)/elasticity-test-transverselyisotropic.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/elasticity-test-transverselyisotropic
	$(COMPILER) $(PATHRHEOTEST)/rheology-test-isotropic-porous.f90        $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHRHEOTEST)/rheology-test-isotropic-porous

aux-tests: $(SPECFAB).o
	$(COMPILER) $(PATHTEST)/reduced-form/reduced-form.f90              $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHTEST)/reduced-form/reduced-form
	$(COMPILER) $(PATHTEST)/rotate-Wigner-D/discrete-rotation-test.f90 $(OBJALL) $(OPTS) $(OPTSNETCDF) -o $(PATHTEST)/rotate-Wigner-D/discrete-rotation-test

#---------------------------

lib$(SPECFAB).so: $(SRCALL) $(OBJALL)
	$(COMPILER) $(OPTS) -shared $(SRCALL) -o $@

$(SPECFAB).o: $(SPECFAB).f90 $(OBJCORE)
	$(COMPILER) $(OPTS) -c $<

$(DEPLESS): %.o: %.f90
	$(COMPILER) $(OPTS) -c $<
	
#---------------------------

$(RHEOLOGIES).o: $(RHEOLOGIES).f90 $(HEADER).o $(TENPROD).o $(LAMBDASOLVER).o
	$(COMPILER) $(OPTS) -c $<

$(ELASTICITIES).o: $(ELASTICITIES).f90 $(HEADER).o $(TENPROD).o
	$(COMPILER) $(OPTS) -c $<

$(DYNAMICS).o: $(DYNAMICS).f90 $(HEADER).o tensorialdynamics.f90 $(TENPROD).o $(MOMENTS).o $(GAUNT).o
	$(COMPILER) $(OPTS) -c $<

$(HOMOGEN).o: $(HOMOGEN).f90 $(HEADER).o $(MOMENTS).o $(RHEOLOGIES).o $(ELASTICITIES).o
	$(COMPILER) $(OPTS) -c $<

$(ENHFAC).o: $(ENHFAC).f90 $(HEADER).o $(MOMENTS).o $(TENPROD).o $(HOMOGEN).o
	$(COMPILER) $(OPTS) -c $<

$(MOMENTS).o: $(MOMENTS).f90 $(HEADER).o $(MANDEL).o
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(FRAMES).o: $(FRAMES).f90 $(MOMENTS).o
	$(COMPILER) $(OPTS) -c $<
	
$(WAVEPROP).o: $(WAVEPROP).f90 $(HEADER).o $(MOMENTS).o $(MANDEL).o $(TENPROD).o $(ELASTICITIES).o $(HOMOGEN).o
	$(COMPILER) $(OPTS) -c $<

$(ROTATION).o: $(ROTATION).f90 $(HEADER).o $(MOMENTS).o 
	$(COMPILER) $(OPTS) -c $<

$(DEFORMMODES).o: $(DEFORMMODES).f90 $(HEADER).o $(TENPROD).o 
	$(COMPILER) $(OPTS) -c $<

$(IDEALSTATE).o: $(IDEALSTATE).f90 $(HEADER).o $(ROTATION).o $(FRAMES).o
	$(COMPILER) $(OPTS) -c $<

$(LAMBDASOLVER).o:
	rm -f $(OBJM77) $(LAMBDASOLVER).o $(LAMBDASOLVER).mod $(PATHM77)/test
	cd $(PATHM77) && $(COMPILER) $(OPTS) -c *.f
	$(COMPILER) $(OPTS) -c $(LAMBDASOLVER).f90
	$(COMPILER) $(OPTS) $(PATHM77)/test.f90 $(LAMBDASOLVER).o $(OBJM77) -o $(PATHM77)/test

$(ELMER).o: $(ELMER).f90
	$(COMPILER) $(OPTS) -c $<

#---------------------------

### Elmer, ...

#elmer: $(ELMERS)

#$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(MODCORE)
#	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@

