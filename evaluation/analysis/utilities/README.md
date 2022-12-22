# utilities

Collection of scripts and instructions for utility scripts.

- compare_two_ncfiles.py: compare two netcdf files and find the maximum difference if the variable name exists in both files and the values are different. This script is useful for diagnosing reproducibility;

- du_rootdir.py: calculate disk usage under a root directory

- count_files_rootdir.py: count files under a root directory

- copy_variable_nc.py: copy variables from one netcdf file to another;

- rm_variables_nc.py: copy one netcdf file to another except the specified variables;

- monthly_average_ncfiles.py: calculate monthly average from hourly or daily data stored in one or multiple netcdf files

- append_netcdf_files.py: catenate variables from multiple files according to time axis and output them to a single file

- find_nearest_point.py: find the nearest grid point index with certain criteria, for example, only land points or certain vegetation category.
