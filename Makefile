# N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

COMPILER=gfortran

LIBNC=-lnetcdff  
LIBNCDIRS=-I/usr/include -L/usr/lib 
FLAGS=-O1 -m64 -mcmodel=medium -ffree-line-length-1024 -Wall -lm -llapack -lblas -lfftw3 $(LIBNCDIRS) $(LIBNC) 

SPECFAB=specfab
MOMENTS=moments
TENPROD=tensorproducts

demo: $(SPECFAB).o $(MOMENTS).o
	$(COMPILER) demo.f90 $(SPECFAB).o $(TENPROD).o $(MOMENTS).o $(FLAGS) -o demo
	@echo "--------------------------\nTo get going, try running (instructions on how to plot the results will follow):"
	@echo "'./demo 3 uc_zz' for uniaxial compression (uc) in the vertical (z) with a nprime=3 grian rheology"
	@echo "'./demo 1 ue_zz' for uniaxial extension (ue) in the vertical (z) with a nprime=1 grian rheology"
	@echo "'./demo 3 ss_xz' for simple shear (ss) along the x--z plane with a nprime=3 grian rheology"
	
$(SPECFAB).o:
	$(COMPILER) -c $(TENPROD).f90 $(FLAGS)
	$(COMPILER) -c $(SPECFAB).f90 $(FLAGS)

$(MOMENTS).o: 
	@echo "*** Compiling structure tensor expressions... this may take some time but only needs to be done once ***"
	$(COMPILER) -c $(MOMENTS).f90 $(FLAGS)

clean:
	rm -f demo $(TENPROD).o $(MOMENTS).o $(SPECFAB).o
	
clear:
	rm -f demo $(SPECFAB).o

