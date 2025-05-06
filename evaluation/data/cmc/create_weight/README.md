1. Create lat/lon for CMC polar grid using test_proj_sphere

this requires some trial and error since the dx/dy is defined in the product at 60N, but this code 
only has either LL or center as the reference. For ~5 decimal point accuracy, use

   proj%dx       = 23813.253175
   proj%dy       = 23813.253175

use center (353,353) as pole for the reference - note that this grid is therefore not symmetric 

   proj%lat1     = 90.0 
   proj%lon1     =  0.0 

also seems the documentation is wrong for this grid, but it may also be the code definition,
set the stdlon = 280.0 (not 10), this could be difference in definition of stdlon

   proj%stdlon   = 280.0

2. create_cmc_scrip.ncl to create a scrip file

3. create a weights file

    - ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C96/fix_20231027/C96_SCRIP.nc \
       --weight CMC_polar_stereo-C96_nearest_wts.nc --method neareststod

    - ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C384/fix_20231027/C384_SCRIP.nc \
       --weight CMC_polar_stereo-C384_nearest_wts.nc --method neareststod

    - ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C768/fix_20231027/C768_SCRIP.nc \
       --weight CMC_polar_stereo-C768_nearest_wts.nc --method neareststod

     - ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C1152/fix_20231027/C1152_SCRIP.nc \
       --weight CMC_polar_stereo-C1152_nearest_wts.nc --method neareststod

sbatch -A fv3-cpu --time=1:00:00 -n 1 --wrap "ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C768/fix_20231027/C768_SCRIP.nc --weight CMC_polar_stereo-C768_nearest_wts.nc --method neareststod --ignore_unmapped"

ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/prototype/prototype_SCRIP.nc \
       --weight CMC_polar_stereo-prototype_nearest_wts.nc --method neareststod --ignore_unmapped

ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/hr/hr_SCRIP.nc \
       --weight CMC_polar_stereo-hr_nearest_wts.nc --method neareststod --ignore_unmapped

sbatch -A fv3-cpu --time=1:00:00 -n 1 --wrap "ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc --destination /scratch2/NCEPDEV/land/data/evaluation/domains/prototype/prototype_SCRIP.nc --weight CMC_polar_stereo-prototype_nearest_wts.nc --method neareststod --ignore_unmapped"

sbatch -A fv3-cpu --time=1:00:00 -n 1 --wrap "ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc --destination /scratch2/NCEPDEV/land/data/evaluation/domains/hr/hr_SCRIP.nc --weight CMC_polar_stereo-hr_nearest_wts.nc --method neareststod --ignore_unmapped"

ESMF_RegridWeightGen --ignore_degenerate --source ./CMC_polar_stereo_SCRIP.nc \
       --destination /scratch2/NCEPDEV/land/data/evaluation/domains/C96_conus/fix_20231027/C96_conus_SCRIP.nc \
       --weight CMC_polar_stereo-C96_conus_nearest_wts.nc --method neareststod --ignore_unmapped

Michael Barlage & Youlong Xia
March 13, 2024
