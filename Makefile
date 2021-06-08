# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2021

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
COMPILER=gfortran -ffree-line-length-none -m64 -Wall -fPIC
OPTS=-O1 -mcmodel=medium -lm -llapack -lblas -lfftw3
OPTSNETCDF=-lnetcdff -I/usr/include -L/usr/lib 

SPECFAB=specfab
MOMENTS=moments
GAUNT=gaunt
TENPROD=tensorproducts

ALLOBJS=$(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(GAUNT).o
ALLSRCS=$(ALLOBJS:.o=.f90)

########################

demo: $(SPECFAB).o
	$(COMPILER) demo/demo_LATROT.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o demo/demo_LATROT
	$(COMPILER) demo/demo_DDRX.f90   $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o demo/demo_DDRX
	mkdir -p demo/solutions
	@echo "\n*** To get going, try the demos:"
	@echo "\n--- Lattice rotation ---"
	@echo "cd demo; ./demo_LATROT uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_LATROT ss_xz :: simple shear (ss) along the x--z plane"
	@echo "\n--- Dynamic recrystallization ---"
	@echo "cd demo; ./demo_DDRX uc_zz :: uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; ./demo_DDRX ss_xz :: simple shear (ss) along the x--z plane"

demoso: lib$(SPECFAB).so
	$(COMPILER) demo/demo_LATROT.f90 -L./ -lspecfab $(OPTS) $(OPTSNETCDF) -o demo/demo_LATROT
	mkdir -p demo/solutions

$(SPECFAB)py: $(SPECFAB).o
	f2py -lm -llapack -lblas -lfftw3 -I. $(ALLOBJS) -c -m specfabpy specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet 
	mkdir -p demo/solutions
	cp specfabpy.cpython* demo/
	@echo "\n*** To get going, try the demo:"
	@echo "cd demo; python3 demo_LATROT.py uc_zz ::: for uniaxial compression (uc) in the vertical (z)"
	@echo "cd demo; python3 demo_LATROT.py ss_xz ::: for simple shear (ss) along the x--z plane"

########################

$(SPECFAB).o: $(MOMENTS).o $(GAUNT).o
	$(COMPILER) $(OPTS) -c $(TENPROD).f90 
	$(COMPILER) $(OPTS) -c $(SPECFAB).f90

lib$(SPECFAB).so: $(ALLSRCS) $(ALLOBJS)
	$(COMPILER) $(OPTS) -shared $(ALLSRCS) -o $@

########################

$(MOMENTS).o: 
	@echo "***************************************************************************************************"
	@echo "*** Compiling structure tensor expressions... this may take some time but is required only once ***"
	@echo "***************************************************************************************************"
	$(COMPILER) $(OPTS) -c $(MOMENTS).f90

$(GAUNT).o: 
	@echo "*****************************************************************************************"
	@echo "*** Compiling gaunt coefficients... this may take some time but is required only once ***"
	@echo "*****************************************************************************************"
	$(COMPILER) $(OPTS) -c $(GAUNT).f90

clean:
	rm -f demo/demo_LATROT demo/demo_DDRX *.o *.mod *.so
	
clear:
	rm -f demo/demo_LATROT demo/demo_DDRX $(SPECFAB).o $(SPECFAB).mod *.so

