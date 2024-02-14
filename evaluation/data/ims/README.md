This is a  directory to regrid 1km and 4km IMS snow cover into various FV3 grid configurations. It includes several sub-directories below:

(1) fv3_mapping is used to produce IMS index files

(1a) Use build.sh to buil executable file
(1b) modify fv3_mapping.nml to put the corresponding FV3 atmos (Cxx) and ocean (mx) combination, such as C96.mx100, as well as various files' pathes so that the executable can find right path and files

For IMS1km,  
source_i_size=24576
source_j_size=24576

for IMS4km
source_i_size=6144
source_j_size=6144

Generally we use otype="C96.mx100_oro_data",  otype="C384.mx025_oro_data",  otype="C768.mx025_oro_data",  otype="C1152.mx025_oro_data" for the four configurations.

(1c) Use "sbatch submit_mapping.sh" to subbit batch job

(2) land_mapping is used to convert index files into land mapping file (1D vector) for snow DA and ufs-land only run

The same steps are used as (1) uses.

(3) land_mapping_openMP is used for calculation of land mapping files when 1km IMS index files are used  for C384 and C768 as land_mapping cannot complete the final run within 8-hour limit

This sub-directory is similar to the land_mapping directory but use OpenMP to compile the code and job submission

(4) gaussian_mapping includes four ncl scripts to process land mapping for prototype and hr (2D gaussian grid)

(5) sorc subdirectory includes two f90 codes to process C96, C384, C768, and C1152 1D vector and prototype and hr 2D gaussian grid. Use "build.sh" to build executables

(6) For C96, C384, C768, C1152, prototype, and hr, there are ims_regrid.nml,  run_ims_snow.sh, and submit_imsRun.sh. In the ims_regrid.nml, IMS4km or IMS1km was set so that the code uses 1km or 4km IMS nc files.

Youlong Xia and Michael Barlage
13 February 2024
