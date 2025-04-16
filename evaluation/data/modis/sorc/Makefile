# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort -132
F90FLAGS       =       -r8 -mcmodel=large 
#F90FLAGS       =
NETCDFMOD      =       -I${NETCDF}/include
NETCDFLIB      =       -L${NETCDF}/lib -lnetcdf -lnetcdff

OBJS =	regrid_modis_cgrid.o regrid_modis_to_gaussian.0

all:	regrid_modis_cgrid.exe regrid_modis_to_gaussian.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(NETCDFMOD) $(*).f90

regrid_modis_cgrid.exe: regrid_modis_cgrid.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_modis_cgrid.o $(NETCDFLIB)

regrid_modis_to_gaussian.exe: regrid_modis_to_gaussian.o
	$(COMPILERF90) -o $(@) $(F90FLAGS) regrid_modis_to_gaussian.o $(NETCDFLIB)

clean:
	rm -f *.exe *.o

