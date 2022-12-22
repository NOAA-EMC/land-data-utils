#!/usr/bin/env python3
#########################################################################################
## Catenate variables from multiple files according to time axis and variables without 
## time axis will be just copied from the first file without any change, examples:
##   1. python append_netcdf_files.py -d source_dir -i input_files -o output_file.nc
##      The input files "input_files*.nc" in the directory "source_dir" will be 
##      appended together into a new file "output_file.nc"
##   2. python append_netcdf_files.py -d source_dir -i input_files_????.nc -o target.nc 
##      -s lon,lat
##      The input files "input_files_????.nc" in the directory "source_dir" will be 
##      appended together into a new file "target.nc", and the variables "lon" and "lat"
##      will be skipped.
##   3. python append_netcdf_files.py -d source_dir -i input_files_????.nc -o target.nc 
##      -v varA,varB
##      Only variables "varA" and "varB" in the input files "input_files_????.nc" in
##      the directory "source_dir" will be copied (if they do not have time axis) or
##      catenated (if they have time axis) into a new file "target.nc". 
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
#########################################################################################
import glob
import argparse
import os
from os.path import exists
import sys
import netCDF4 as nc
import numpy as np

def last_24chars(x):
  return(x[-14:])

def append_ncfiles(rootd, ifname, ofname, vname, skip):
  skiplist = skip.split(",")
  varlist = vname.split(",")
  if not rootd == '':
    if not ifname.endswith('.nc'):
      all_file  = glob.glob(rootd + '/' + ifname + '*.nc')
    else:
      all_file  = glob.glob(rootd + '/' + ifname)
  else:
    if not ifname.endswith('.nc'):
      all_file  = glob.glob(ifname + '*.nc')
    else:
      all_file  = glob.glob(ifname)
  all_files = sorted(all_file, key=last_24chars)
  print("-------------------------")
  print("The files to be appended:")
  print("    First: ", all_files[0])
  print("    Last:  ", all_files[len(all_files)-1])
# for file in all_files:
#   print("    ", file)
  ifname = all_files[0]
  with nc.Dataset(ifname) as src:
    with nc.Dataset(ofname, 'w', format=src.file_format) as dst:
      # copy global attributes all at once via dictionary
      dst.setncatts(src.__dict__)
      # copy dimensions
      for name, dimension in src.dimensions.items():
        dst.createDimension(
          name, (len(dimension) if not dimension.isunlimited() else None))
      # copy all file data except for the excluded
      vnum = 0
      for name, variable in src.variables.items():
        if name in skiplist:
          continue
        if not vname == '' and not name in varlist:
          continue
        vnum += 1
        print("Variable: ", vnum, name)
        createattrs = variable.filters()
        if createattrs is None:
          createattrs = {}
        else:
          chunksizes = variable.chunking()
          if chunksizes == "contiguous":
            createattrs["contiguous"] = True
            print("contiguous: ",createattrs["contiguous"])
          else:
            createattrs["chunksizes"] =  chunksizes
        x = dst.createVariable(name, variable.datatype, variable.dimensions, **createattrs)
        # copy variable attributes all at once via dictionary
        dst[name].setncatts(src[name].__dict__)
        if 'time' in variable.dimensions:
          num = 0
          for fname in all_files:
            infile = nc.Dataset(fname, "r")
            num += 1
#           print(num, name, fname)
            data = infile.variables[name][...]
            if num == 1:
              var = data
            else:
              var = np.append(var, data, axis=0)
            infile.close()
          dst[name][:] = var
        else:
          dst[name][:] = src[name][:]
  print("The output file: ", ofname)
  print("The script ended normally!")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--rootd',  help="name of the root directory", default="")
    ap.add_argument('-i', '--ifname', help="input files or their prefix", required=True)
    ap.add_argument('-o', '--ofname', help="output file name", required=True)
    ap.add_argument('-v', '--vname',  help="variables to be appended", default="")
    ap.add_argument('-s', '--skip',   help="skip list", default="")
    MyArgs = ap.parse_args()
    append_ncfiles(MyArgs.rootd, MyArgs.ifname, MyArgs.ofname, MyArgs.vname, MyArgs.skip)
