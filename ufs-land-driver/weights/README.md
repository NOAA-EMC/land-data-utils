# weights

Modify `create_weights.sh` for your case. Submit script using:

`sbatch create_weights.sh`

Some standard grids may already be created here:

`/scratch4/NCEPDEV/land/data/ufs-land-driver/weights`

 set parameters for weights generation
 
`atm_res`      : global or regional grid resolution [C48,C96,C192,C384,C768,C1152] 

`ocn_res`      : ocean resolution [mx500,mx100,mx050,mx025], not used for regional

`grid_version` : fix file version [20231027 (global grids consistent with GFSv17), AQM (regional air quality grid), ARC (regional Arctic grid)] 

`fixfile_path` : top level path for fix files, e.g., `/scratch3/NCEPDEV/global/role.glopara/fix/orog/` 

`grid_extent`  : grid options [total,subset] 
 
`subset_name`  : if a subset grid, name for subset, e.g., conus

`datm_source`  : data atmosphere source [ERA5,GDAS,CDAS]

`datm_source_file` : a single datm_source file, e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/2022/ERA5_forcing_2022-12-31.nc`

`destination_scrip_path` : destination SCRIP path, e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/`

Outputs in current directory `atm_res`.`ocn_res` (`C96.mx100` example):

`ERA5-C96.mx100_bilinear_wts.nc` :
* contains ESMF weights to regrid gridded ERA5 to C96.mx100 vector using bilinear interpolation

`ERA5-C96.mx100_neareststod_wts.nc` :
* contains ESMF weights to regrid gridded ERA5 to C96.mx100 vector using bilinear interpolation

