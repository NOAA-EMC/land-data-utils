#

OBJS = cmc_proj_sphere.o ps_routines.o
 NETCDFMOD      =       -I${NETCDF}/include
 NETCDFLIB      =       -L${NETCDF}/lib -lnetcdf -lnetcdff

all : $(OBJS)
	ifort -r8 -o cmc_proj $(NETCDFMOD)  $(OBJS) $(NETCDFLIB)


cmc_proj_sphere.o : cmc_proj_sphere.f90
	ifort -r8 $(NETCDFMOD) -c cmc_proj_sphere.f90


ps_routines.o : ps_routines.f90
	ifort -r8 $(NETCDFMOD) -c ps_routines.f90


cmc_proj_sphere.o: ps_routines.o


clean :
	rm -f *.o *.mod cmc_proj
