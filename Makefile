# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

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

demopy: $(SPECFAB)py
	mkdir -p solutions
	@echo "-----------------------------------------------"
	@echo "To get going, try running (instructions on how to plot the results will follow):"
	@echo "python3 demo.py 20 3 uc_zz ::: for uniaxial compression (uc) in the vertical (z) with a nprime=3 grian rheology and truncation L=20"
	@echo "python3 demo.py 20 1 ue_zz ::: for uniaxial extension (ue) in the vertical (z) with a nprime=1 grian rheology and truncation L=20"
	@echo "python3 demo.py 20 3 ss_xz ::: for simple shear (ss) along the x--z plane with a nprime=3 grian rheology and truncation L=20"

demo: $(SPECFAB).o
	$(COMPILER) demo.f90 $(ALLOBJS) $(OPTS) $(OPTSNETCDF) -o demo
	mkdir -p solutions
	@echo "-----------------------------------------------"
	@echo "To get going, try running (instructions on how to plot the results will follow):"
	@echo "./demo 20 3 uc_zz ::: for uniaxial compression (uc) in the vertical (z) with a nprime=3 grian rheology and truncation L=20"
	@echo "./demo 20 1 ue_zz ::: for uniaxial extension (ue) in the vertical (z) with a nprime=1 grian rheology and truncation L=20"
	@echo "./demo 20 3 ss_xz ::: for simple shear (ss) along the x--z plane with a nprime=3 grian rheology and truncation L=20"

demoso: lib$(SPECFAB).so
	$(COMPILER) demo.f90 -L./ -lspecfab $(OPTS) $(OPTSNETCDF) -o demo
	mkdir -p solutions

########################

$(SPECFAB)py: $(SPECFAB).o
	f2py -lm -llapack -lblas -lfftw3 -I. $(ALLOBJS) -c -m specfabpy specfabpy.f90 --f90flags="-ffree-line-length-none -mcmodel=medium" --quiet
	
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
	rm -f demo *.o *.mod *.so
	
clear:
	rm -f demo $(SPECFAB).o $(SPECFAB).mod *.so

