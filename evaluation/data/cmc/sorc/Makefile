# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

 COMPILERF90    =       ifort -132
 FREESOURCE     =       #-ffree-form  -ffree-line-length-none
 F90FLAGS       =       -r8
 NETCDFMOD      =       -I/apps/netcdf/4.7.0/intel/18.0.5.274/include
 NETCDFLIB      =       -L/apps/netcdf/4.7.0/intel/18.0.5.274/lib -lnetcdf -lnetcdff


OBJS =	regrid_CMCsnow_to1Dvector.o regrid_CMCsnow_toGaussianGrid.o 

all:	regrid_CMCsnow_to1Dvector.exe regrid_CMCsnow_toGaussianGrid.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

regrid_CMCsnow_to1Dvector.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_CMCsnow_to1Dvector.o $(NETCDFLIB)

regrid_CMCsnow_toGaussianGrid.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_CMCsnow_toGaussianGrid.o $(NETCDFLIB)

clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#
