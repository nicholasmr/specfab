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

DEMOPATH=../demo

ALLOBJS=$(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(GAUNT).o
ALLMODS=$(ALLOBJS:.o=.mod)
ALLSRCS=$(ALLOBJS:.o=.f90)

#-------------------
# Elmer, JOSEF, ...
#-------------------

ELMERSRCS=$(wildcard $(DEMOPATH)/elmer*.f90)
ELMERS=$(ELMERSRCS:.f90=)

########################################

#-------------------
# Core
#-------------------

demo: $(SPECFAB).o
	$(COMPILER) $(DEMOPATH)/demo_LATROT.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(DEMOPATH)/demo_LATROT
	$(COMPILER) $(DEMOPATH)/demo_DDRX.f90   $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(DEMOPATH)/demo_DDRX
	mkdir -p $(DEMOPATH)/solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd $(DEMOPATH); ./demo_LATROT uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(DEMOPATH); ./demo_LATROT ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd $(DEMOPATH); ./demo_DDRX uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(DEMOPATH); ./demo_DDRX ss_xz :: simple shear (ss) along the x--z plane"

specfabpy: $(SPECFAB).o
	f2py -lm -llapack -lblas -lfftw3 -I. $(ALLOBJS) -c -m specfabpy specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet 
	mkdir -p $(DEMOPATH)/solutions
	cp specfabpy.cpython* $(DEMOPATH)/
	@echo "\n*** To get going, try the demo:"
	@echo "cd $(DEMOPATH); python3 demo_LATROT.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd $(DEMOPATH); python3 demo_LATROT.py ss_xz ::: for simple shear (ss) along the x--z plane"

lib$(SPECFAB).so: $(ALLSRCS) $(ALLOBJS)
	$(COMPILER) $(OPTS) -shared $(ALLSRCS) -o $@

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

$(MOMENTS).o: 
	@echo "*** Compiling structure tensor expressions and gaunt coefficients... this may take some time but is required only once ***"
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(GAUNT).o: 
	$(COMPILER) $(OPTS) -c $(GAUNT).f90

clean:
	rm -f $(DEMOPATH)/demo_LATROT $(DEMOPATH)/demo_DDRX *.o *.mod *.so
	
clear:
	rm -f $(DEMOPATH)/demo_LATROT $(DEMOPATH)/demo_DDRX $(SPECFAB).o $(SPECFAB).mod *.so

rheologytests: $(SPECFAB).o
	$(COMPILER) $(DEMOPATH)/test_orthotropic_rheology.f90   $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(DEMOPATH)/test_orthotropic_rheology
	$(COMPILER) $(DEMOPATH)/test_tranisotropic_rheology.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $(DEMOPATH)/test_tranisotropic_rheology
	
#-------------------
# Elmer, JOSEF, ...
#-------------------

#elmer: $(ELMERS)

#$(DEMO)/elmer%: $(DEMO)/elmer%.f90 | lib$(SPECFAB).so $(ALLMODS)
#	$(COMPILER) $< $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o $@


