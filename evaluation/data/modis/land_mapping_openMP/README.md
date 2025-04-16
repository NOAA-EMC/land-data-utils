Code uses IMS FV3-mapping index files to create land mapping files for processing 1D vector regrid output (for use of snow DA and ufs-land only run). Input is IMS FV3-mapping index file and output is land maaping file. This is openMP F90 code for 0.05 degree modis and C1152 as land_mapping cannot finish the run task within 8-hours limit

(1) To compile, do "make"

(2) To run, submit executable (in big-mem queue) with submit_mapping.sh, do:
sbatch submit_mapping.sh

Michael Barlage and Youlong Xia.
16 February 2024  
