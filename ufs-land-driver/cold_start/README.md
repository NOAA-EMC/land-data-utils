# cold start initial conditions

To conduct a simulation from a somewhat arbitary set of initial conditions, run this script to get the cold start initial condition file based on:
- grid definition
- source data for temperature
- simulation start date

The date stamp on the file will be one timestep before the simulation start date.

Modify `create_cold_start.sh` for your case. Submit script using:

`sbatch create_cold_start.sh`

Some standard grids and times may already be created here:

`/scratch4/NCEPDEV/land/data/ufs-land-driver/cold_start`

 set parameters for weights generation
 
`atm_res`      : global or regional grid resolution [C48,C96,C192,C384,C768,C1152] 

`ocn_res`      : ocean resolution [mx500,mx100,mx050,mx025], not used for regional

`grid_version` : fix file version [20231027 (global grids consistent with GFSv17), AQM (regional air quality grid), ARC (regional Arctic grid)] 

`grid_extent`  : grid options [total,subset] 
 
`subset_name`  : if a subset grid, name for subset, e.g., conus

`datm_source`  : data atmosphere source [ERA5,GDAS,CDAS]

`datm_source_path` : base path where datm already exists, e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/` 

`static_file_path` : ufs-land-driver static file path (created in Step 1), e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/` 

`yyyy` : year of desired cold start initial condition, e.g., 1999 

`mm` : month of desired cold start initial condition, e.g., 1 

`dd` : day of desired cold start initial condition, e.g., 1 

`hh` : hour of desired cold start initial condition, e.g., 0 

`timestep` : model timestep in minutes, e.g., 60 

Outputs in current directory `atm_res`.`ocn_res` (`C96.mx100` example):

`ERA5-C96.mx100_cold_start_1979-12-31_23:00:00.nc`:
* Cold start IC from ERA5 regrid to C96.mx100 vector for simulation to start 1980-01-01 00:00:00
