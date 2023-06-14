#!/usr/bin/env python3
#########################################################################################
## Calculate monthly average from hourly or daily data stored in one or multiple netcdf 
## files, examples:
##   1. python monthly_average_ncfiles.py -i input_files -o output_file.nc
##      The monthly average will be calculated at each grid for all time-series variables
##      in the input files "input_files*.nc", and output to the file "output_file.nc"
##   2. python monthly_average_ncfiles.py -i input_file.nc -o target.nc -s lon,timestamp
##      The monthly average will be calculated at each grid for data stored in the input
##      file "input_file.nc" for all variables except "lon" and "timestamp", and output
##      to the file "target.nc"
##   3. python monthly_average_ncfiles.py -d dir -i input_files -o target.nc -v sm,tair
##      The monthly average will be calculated at each grid for data stored in the input 
##      files "input_files*.nc" in the directory "dir" for variables "sm" and "tair", 
##      and output to the file "target.nc"
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
#########################################################################################
import glob
import argparse
import os
import sys
import netCDF4 as nc
import numpy as np
from datetime import timedelta
from datetime import datetime

def cal_average(vname, all_files):
  num = 0
  for fname in all_files:
    infile = nc.Dataset(fname, "r")
    num += 1
    data = infile.variables[vname][...]
    if num == 1:
      var = data
    else:
      var = np.append(var, data, axis=0)
    infile.close()
  ovar = np.nanmean(var,axis=0)
  return ovar

def last_14chars(x):
  return(x[-14:])

def convert_timeStr2timeObj(date_time_str):
# formt = '%d/%m/%y %H:%M:%S' # '18/09/19 01:55:19'
  try:
    formt = '%Y-%m-%d %H:%M:%S' # '2016-01-01 00:50:00'
    date_time_obj = datetime.strptime(date_time_str, formt)
  except:
    formt = '%Y-%m-%d %H:%M' # '2016-01-01 00:50:00'
    date_time_obj = datetime.strptime(date_time_str, formt)
  return date_time_obj

def advance_dateTime(date_time_obj, units, dt):
  if units.upper() == 'HOURS':
    return date_time_obj + timedelta(seconds=float(dt*3600))
  elif units.upper() == 'SECONDS':
    return date_time_obj + timedelta(seconds=dt)
  elif units.upper() == 'DAYS':
    return date_time_obj + timedelta(days=int(dt))
  else:
    print("invalid time units")

def average_ncfiles(rootd, ifname, ofname, vname, skip):
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
  all_files = sorted(all_file, key=last_14chars)
  print("-------------------------")
  print("The files to be averaged:")
  if all_files[0] == all_files[len(all_files)-1]:
    print("    File:  ", all_files[0])
  else:
    print("    First: ", all_files[0])
    print("    Last:  ", all_files[len(all_files)-1])

  infile_first = nc.Dataset(all_files[0], "r")
  infile_last  = nc.Dataset(all_files[len(all_files)-1], "r")
  units0 = infile_first.variables['time'].getncattr('units')
  units1 = infile_last.variables['time'].getncattr('units')
  if not units0 == units1:
    sys.exit('Error: different time units')
  time0  = infile_first.variables['time'][:]
  time1  = infile_last.variables['time'][:]
  infile_first.close()
  infile_last.close()
  words   = units0.split(' ')
  time_units = words[0]
  time_begin = words[2]+' '+words[3]
  time_ini = convert_timeStr2timeObj(time_begin)
  time_beg = advance_dateTime(time_ini, time_units, time0[0])
  time_end = advance_dateTime(time_ini, time_units, time1[len(time1)-1])
  year0  = time_beg.year
  month0 = time_beg.month
  year1  = time_end.year
  month1 = time_end.month
  months = (year1-year0)*12 + month1 - month0 + 1
  time_cur = convert_timeStr2timeObj(str(year0)+'-'+str(month0)+'-15 00:00:00')
  difference = time_cur - time_ini
  if months > 1:
    sys.exit("Error: the files contain data more than one month!")

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
      for name, variable in src.variables.items():
        if name in skiplist:
          continue
        if not vname == '' and  not name in varlist:
          continue
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
        if name == 'time':
          dst[name].long_name = 'time'
          dst[name].units = 'days since ' + words[2] + ' ' + words[3]
        else:
          dst[name].setncatts(src[name].__dict__)
        if 'time' in variable.dimensions:
          if name == 'time':
            dst[name][0] = difference.days
          else:
            num = 0
            for fname in all_files:
              infile = nc.Dataset(fname, "r")
              num += 1
              data = infile.variables[name][...]
              if num == 1:
                var = data
              else:
                var = np.append(var, data, axis=0)
              infile.close()
            ovar = np.nanmean(var,axis=0)  
            dst[name][0] = ovar
        else:
          dst[name][:] = src[name][:]
  print("The output file: ", ofname)
    
if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-i', '--ifname', help="input files or their prefix", required=True)
  ap.add_argument('-o', '--ofname', help="output file name", required=True)
  ap.add_argument('-d', '--rootd',  help="name of the root directory", default="")
  ap.add_argument('-v', '--vname',  help="variables to be averaged", default="")
  ap.add_argument('-s', '--skip',   help="skip list", default="timestep")
  MyArgs = ap.parse_args()
  average_ncfiles(MyArgs.rootd, MyArgs.ifname, MyArgs.ofname, MyArgs.vname, MyArgs.skip)
