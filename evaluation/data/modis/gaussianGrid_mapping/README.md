This ncl to produce gaussian grid mapping. As modis lon is from -180 to 180 and corners lon is from 0 to 360, ncl script has shifted to fit 0 - 360 degree. The FORTRAN 90 code can read mapping file and modis file directly.

The ncl scripts are in evaluation/data/modis/gaussianGrid_mapping/. To run the ncl scripts, you need to load ncl module first:

module load ncl
ncl create_pt_mapping_hr.ncl
ncl create_pt_mapping_prototype.ncl


Youlong Xia and Michael Barlage
