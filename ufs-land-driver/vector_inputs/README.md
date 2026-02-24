# vector_inputs

Modify `create_inputs.sh` for your case. Submit script using:

`sbatch create_inputs.sh`

Some standard grids may already be created here:

`/scratch4/NCEPDEV/land/data/ufs-land-driver/vector_inputs`

 set parameters for grid generation
 
`atm_res`      : global or regional grid resolution [C48,C96,C192,C384,C768,C1152] 

`ocn_res`      : ocean resolution [mx500,mx100,mx050,mx025], not used for regional

`grid_version` : fix file version [20231027 (global grids consistent with GFSv17), AQM (regional air quality grid), ARC (regional Arctic grid)] 

`fixfile_path` : top level path for fix files, e.g., `/scratch3/NCEPDEV/global/role.glopara/fix/orog/` 

`grid_extent`  : grid options [total,subset] 
 
`subset_name`  : if a subset grid, name for subset, e.g., conus

`subset_maxlat`  : cutout maximum latitude

`subset_minlat`  : cutout minimum latitude

`subset_maxlon`  : cutout maximum longitude

`subset_minlon`  : cutout minimum longitude

Outputs in current directory (`C96.mx100` example):

`ufs-land_C96.mx100_corners.nc` : 
* contains lat/lon of grid centers and corners organized in vector
* vector is organized from tile 1 to 6 starting in lower-left corner with x-dimension faster-varying
* check file metadata for grid information

`ufs-land_C96.mx100_static_fields.nc` :
* contains fix file inputs put into vector format, a necessary input for the ufs-land-driver

`ufs-land_C96.mx100_SCRIP.nc` :
* contains SCRIP "unstructured" format information used for ESMF regridding

`ufs-land_C96.mx100_SCRIP_veg.nc` :
* contains SCRIP "unstructured" format information for vegetated grids used for ESMF regridding
* this and the following two are used for regridding vectors between resolutions to make certain grids are consistently matched (veg->veg, bare->bare, snow/ice->snow/ice)

`ufs-land_C96.mx100_SCRIP_bare.nc`   :
* contains SCRIP "unstructured" format information for bare grids used for ESMF regridding

`ufs-land_C96.mx100_SCRIP_snow.nc`   :
* contains SCRIP "unstructured" format information for snow/ice grids used for ESMF regridding
