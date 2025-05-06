CMC snow depth

The CMC ascii snow depth data can by downloaded using the following scripts: https://github.com/nsidc/nsidc0447-scripts

Daily snow depth for 2012-2020

Data are on ~2km stereographic grid

To regrid to different grids, use the following steps:

1. see https://github.com/NOAA-EMC/land-data-utils/tree/main/ufs-land-driver/weights to create weights. 
Only needs to be done once when fv3 grid and interpolation method are determined.

2. On Hera, need to load modules first:
module use /scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/modulefiles/Core

module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load netcdf/4.7.0
module load netcdf-hdf5parallel/4.7.4

3. Enter sorc directory and do "make" to compiled f90 code to produce excutables

4. Enter each individual destination grid directory such as C96_conus, C96, C384, C768, C1152, prototype, and hr, do "sbatch submit_cmcRun.sh"

Noted that destination_locs and weight_locs are directly read from weights files, such as /scratch2/NCEPDEV/land/data/evaluation/CMC/fix_20231027/weigh
ts/CMC_polar_stereo-C96_nearest_wts.nc, using ncdump -h to check this file to find the corresponding values below.

destination_locs=n_b 
weight_locs=n_s 

These values depend on fv3 grid and interpolation method. 

Youlong Xia and Michael Barlage
July 31, 2024
