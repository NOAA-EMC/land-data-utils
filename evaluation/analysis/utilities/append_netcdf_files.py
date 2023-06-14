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
##   4. python append_netcdf_files.py -d source_dir -i input_files -o output_file.nc
##      -t 2011100100-2012093023
##      The input files "input_files*.nc" in the directory "source_dir" will be 
##      appended together into a new file "output_file.nc" only when data records
##      during period 2011-10-01 00:00 and 2012-09-30 23:00 are found in those files.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
#########################################################################################
import glob
import argparse
import os
from os.path import exists
import sys
import netCDF4 as nc
import numpy as np
from datetime import timedelta
from datetime import datetime

def last_24chars(x):
  return(x[-24:])

def append_ncfiles(rootd, ifname, ofname, vname, skip, trange):
  skiplist = skip.split(",")
  varlist = vname.split(",")
  if not trange == '':
    tranges = trange.split('-')
    time_beg = datetime.strptime(tranges[0], '%Y%m%d%H')
    time_end = datetime.strptime(tranges[1], '%Y%m%d%H')
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
  flag = np.ones(len(all_files), dtype=np.int)
  first_file = None
  last_file  = None
  if not trange == '':
    for fid in range(len(all_files)):
      fname = all_files[fid]
      infile     = nc.Dataset(fname, "r")
      times      = infile.variables['time'][:]
      units      = infile.variables['time'].getncattr('units')
      words      = units.split(' ')
      time_units = words[0]
      time_ini   = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M:%S')
      factor = 1
      if time_units.upper() == 'DAYS':
        factor = 86400
      elif time_units.upper() == 'HOURS':
        factor = 3600
      elif not time_units.upper() == 'SECONDS':
        sys.exit("Error: time units not recognized!")
      time_start = time_ini + timedelta(seconds=times[0]*factor)
      time_stop  = time_ini + timedelta(seconds=times[len(times)-1]*factor)
      if time_start > time_end or time_stop < time_beg:
        flag[fid] = 0
      if flag[fid] > 0:
        if first_file == None:
          first_file = fname
        last_file = fname
      infile.close()
  else:
    first_file = all_files[0]
    last_file  = all_files[len(all_files)-1]
  print("The files to be appended:")
  print("    First: ", first_file)
  print("    Last:  ", last_file)
# for file in all_files:
#   print("    ", file)
  ifname = all_files[0]
  with nc.Dataset(ifname) as src:
    with nc.Dataset(ofname, 'w', format=src.file_format) as dst:
      # copy global attributes all at once via dictionary
      dst.setncatts(src.__dict__)
      # copy dimensions
      for name, dimension in src.dimensions.items():
        if name.upper() == 'TIME':
          dst.createDimension(name, None)
        else:
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
        try:
          the_attrs = src[name].__dict__
          aflag = 0
          afill = None
          for aname, avalue in the_attrs.items():
            if aname == '_FillValue':
              aflag = 1
              afill = avalue
          if aflag == 0:
            dst[name].setncatts(src[name].__dict__)
          else:
            print(name, 'Fille: ', afill)
#           dst[name]._FillValue = afill
            dst[name].FillValue = afill
            print(name, 'Fille: ', afill,' done')
        except:
          print("Source atts for ",name,": ",src[name].__dict__)
        if 'time' in variable.dimensions:
          num = 0
          for fid in range(len(all_files)):
            if flag[fid] > 0:
              fname = all_files[fid]
              infile = nc.Dataset(fname, "r")
              num += 1
#             print(num, name, fname)
              data = infile.variables[name][...]
              if num == 1:
                var = data
              else:
                var = np.append(var, data, axis=0)
              infile.close()
          if num < 1:
            sys.exit("Error: no data record found!")
          dst[name][:] = var
        else:
          dst[name][:] = src[name][:]
  print("The output file: ", ofname)
  print("The script ended normally!")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--rootd',  help="name of the root directory",  default="")
    ap.add_argument('-i', '--ifname', help="input files or their prefix", required=True)
    ap.add_argument('-o', '--ofname', help="output file name",            required=True)
    ap.add_argument('-v', '--vname',  help="variables to be appended",    default="")
    ap.add_argument('-s', '--skip',   help="skip list",                   default="")
    ap.add_argument('-t', '--trange', help="time range for plotting",     default="")
    MyArgs = ap.parse_args()
    append_ncfiles(MyArgs.rootd, MyArgs.ifname, MyArgs.ofname, MyArgs.vname, MyArgs.skip, MyArgs.trange)
