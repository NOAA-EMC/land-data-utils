# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

COMPILERF90    =       ifort
NETCDFMOD      =       -I${netcdf_fortran_ROOT}/include
NETCDFLIB      =       -L${netcdf_fortran_ROOT}/lib -lnetcdff

all:	create_core_radiation_from_monthly.exe create_shell_daily_files.exe create_core_hourly_from_monthly.exe 

create_core_hourly_from_monthly.exe: create_core_hourly_from_monthly.f90
	$(COMPILERF90) -o $(@) -diag-disable=10448 \
	 $(NETCDFMOD) create_core_hourly_from_monthly.f90 $(NETCDFLIB)

create_shell_daily_files.exe: create_shell_daily_files.f90
	$(COMPILERF90) -o $(@) -diag-disable=10448 \
	 $(NETCDFMOD) create_shell_daily_files.f90 $(NETCDFLIB)

create_core_radiation_from_monthly.exe: create_core_radiation_from_monthly.f90
	$(COMPILERF90) -o $(@) -diag-disable=10448 \
	 $(NETCDFMOD) create_core_radiation_from_monthly.f90 $(NETCDFLIB)

clean:
	rm -f *.exe *.o *.mod
