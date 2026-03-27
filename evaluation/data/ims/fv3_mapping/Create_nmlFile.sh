#! /bin/sh
# IMS data: 1km - 24576, 4km - 6144, 24km - 1024
# IMS data source: IMS.4km or IMS.1km
# FV3 grid resolution and otype: 96 -> C96.mx100_oro_data, 
# 384->C384.mx025_oro_data, 768->C768.mx025_oro_data
# 

fv3Resolution=96
imsSize=24576
imsResolution=1km
orogDirectory=/scratch1/NCEPDEV/global/glopara/fix/orog/20231027
tileDirectory=/scratch1/NCEPDEV/global/glopara/fix/orog/20231027
imscordDirectory=/scratch2/NCEPDEV/land/data/evaluation/IMS/fix_coords/
outputDirectory=/scratch2/NCEPDEV/land/data/evaluation/IMS/fix_20231027/index_files/

cat > fv3_mapping.nml << EOF
&fv3_mapping_nml
 tile_dim=$fv3Resolution
 otype="C${fv3Resolution}.mx100_oro_data"
 orog_path="${orogDirectory}/C${fv3Resolution}"
 tile_path="${tileDirectory}/C${fv3Resolution}"
 obs_source="IMS${imsResolution}"
 source_i_size=$imsSize
 source_j_size=$imsSize
 ims_path="${imscordDirectory}"
 ims_lat_name="imslat_${imsResolution}_8bytes.bin"
 ims_lon_name="imslon_${imsResolution}_8bytes.bin"
 out_path="${outputDirectory}"
/
EOF

