# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

#COMPILER=ifort -free -xHost -shared-intel -align all  -debug all -qopt-report-phase=vec -diag-disable 8291 -diag-disable 8290  
#INCLUDES=-I/usr/local/include -L/usr/local/lib

COMPILER=gfortran -ffree-line-length-1024 -m64 -Wall
INCLUDES=-I/usr/include -L/usr/lib 

OPTS=-O1 -mcmodel=medium -lm -llapack -lblas -lfftw3 -lnetcdff $(INCLUDES)

SPECFAB=specfab
MOMENTS=moments
TENPROD=tensorproducts

demo: $(MOMENTS).o $(SPECFAB).o
	$(COMPILER) demo.f90 $(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(OPTS) -o demo
	@echo "-----------------------------------------------"
	@echo "To get going, try running (instructions on how to plot the results will follow):"
	@echo "'./demo 3 uc_zz' for uniaxial compression (uc) in the vertical (z) with a nprime=3 grian rheology"
	@echo "'./demo 1 ue_zz' for uniaxial extension (ue) in the vertical (z) with a nprime=1 grian rheology"
	@echo "'./demo 3 ss_xz' for simple shear (ss) along the x--z plane with a nprime=3 grian rheology"
	
$(SPECFAB).o:
	$(COMPILER) -c $(TENPROD).f90 $(OPTS)
	$(COMPILER) -c $(SPECFAB).f90 $(OPTS)

$(MOMENTS).o: 
	@echo "*** Compiling structure tensor expressions... this may take some time but only needs to be done once ***"
	$(COMPILER) -c $(MOMENTS).f90 $(OPTS)

clean:
	rm -f demo *.o *.mod
	
clear:
	rm -f demo $(SPECFAB).o $(SPECFAB).mod

