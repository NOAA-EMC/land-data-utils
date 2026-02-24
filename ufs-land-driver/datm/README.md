# create data atmosphere (datm) files

Modify `create_datm_files.sh` for your case. Submit script using:

`sbatch create_datm_files.sh`

Some standard grids may already be created here:

`/scratch4/NCEPDEV/land/data/ufs-land-driver/datm`

 set parameters for datm generation
 
`atm_res`      : global or regional grid resolution [C48,C96,C192,C384,C768,C1152] 

`ocn_res`      : ocean resolution [mx500,mx100,mx050,mx025], not used for regional

`grid_version` : fix file version [20231027 (global grids consistent with GFSv17), AQM (regional air quality grid), ARC (regional Arctic grid)] 

`grid_extent`  : grid options [total,subset] 
 
`subset_name`  : if a subset grid, name for subset, e.g., conus

`datm_source`  : data atmosphere source [ERA5,GDAS,CDAS,CORe]

`datm_source_path` : datm_source files path, e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/ `

`elevation_source_filename` : source data elevation file, e.g., `/scratch4/NCEPDEV/land/data/ufs-land-driver/datm/ERA5/original/elevation/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc`

`static_file_path` : ufs-land-driver static file path (created in Step 1), e.g., `"/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs/"` 

`weights_path` : weights files path (created in Step 2), e.g., `"/scratch4/NCEPDEV/land/data/ufs-land-driver/weights/"`

`interpolation_method1` : primary ESMF regrid method [bilinear,neareststod]; used for all inputs except optionally precipitation

`interpolation_method2` : secondary ESMF regrid method [bilinear,neareststod]; optional method used for precipitation

`precip_interpolation_method` : ESMF regrid method [all,method1,method2]

`regrid_tasks_file` : file with all the individual commands to create datm files, e.g., `"regrid-tasks.ERA5"` 

Outputs in current directory `atm_res`.`ocn_res` (`C96.mx100` example):

`ERA5-C96.mx100_datm_2009-07-01.nc` :
* atmospheric forcing from ERA5 regrid to C96.mx100 vector for 2009-07-01
