Code to create file with indexes used to map from ASCII IMS file onto FV3 grid. 

To compile on hera, use build.sh. 

requires following modules: 

module use /scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/u
nified-env-rocky8/install/modulefiles/Core

module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load netcdf/4.7.0
module load netcdf-hdf5parallel/4.7.4

To run, submit executable (in big-mem queue)

There is a testcase and example submission script for C48 on hera at: 
/scratch2/BMC/gsienkf/Clara.Draper/DA_test_cases/snow/IMSobsproc/get_index/

This code is slightly modified to create mapping file for MODIS data.

The potential FV3 atmos and ocean configurations are: C48mx500, C48mx100, C96mx500, C96mx100, C192mx025, C192mx100, C384mx025, C768mx025, C1152mx025 

For offline snow DA workflow, we usually use C48mx100, C96mx100, C384mx025, C768mx025, C1152mx025

Mike Barlage, Clara Draper, Youlong Xia.
