CMC snow depth

The CMC ascii snow depth data can by downloaded using the following scripts: https://github.com/nsidc/nsidc0447-scripts

Daily snow depth for 2012-2020

Data are on ~2km stereographic grid

To regrid to different grids, use the following steps:

1. enter to /scratch2/NCEPDEV/land/Youlong.Xia/land-data-utils/evaluation/data/cmc/create_weight to produce weights data. See the README.md in that directory. 
Only needs to be done once when fv3 grid and interpolation method are determined.

2. Enter sorc directory and do "make" to compiled f90 code to produce excutables

3. Enter each individual fv3 grid such as C96_conus, C96, C384, C768, C1152, prototype, and hr, do "sbatch submit_cmcRun.sh"

Noted that destination_locs and weight_locs are directly read from weights files, such as /scratch2/NCEPDEV/land/data/evaluation/CMC/fix_20231027/weigh
ts/CMC_polar_stereo-C96_nearest_wts.nc, using ncdump -h to check this file to find the corresponding values below.

destination_locs=n_b 
weight_locs=n_s 

These values depend on fv3 grid and interpolation method. 

Youlong Xia and Michael Barlage
March 14, 2024
