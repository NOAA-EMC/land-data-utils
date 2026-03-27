The ncl scripts are used to create prototype and hr mapping files for use of regridding ims snow cover fraction into gaussian grid.

IMS4km_prototype_mapping.nc, IMS4km_hr_mapping.nc, IMS4km_hr_mapping.nc, IMS1km_hr_mapping.nc

prototype is for C384 Gaussian grid and hr is for C768 Gaussian grid. These mapping files are used to cound how many ims grid within a given FV3 grid and calculate their average value in F90 code.

The scripts can be run in interactive mode or bash script. However, create_pt_ims1km_mapping_hr.ncl cannot complete its run in 8 hours maximum limit. It is needed to run create_pt_ims1km_mapping_hr.ncl interactively.

Youlong Xia and Michael Barlage

July 31, 2024
