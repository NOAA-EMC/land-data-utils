#!/usr/bin/env python3
#########################################################################################
## compare two netcdf files and find the maximum difference if the variable name exists
## in both files and the values are different. This script is useful for diagnosing 
## reproducibility. Examples:
##   1. python compare_two_ncfiles.py -if input_fileA.nc,input_fileB.nc
##      Compare two netcdf files "input_fileA.nc" and "input_fileB.nc"
##   2. python compare_two_ncfiles.py -rd dir -if input_fileA.nc,input_fileB.nc
##      Compare two netcdf files "input_fileA.nc" and "input_fileB.nc" which are located
##      under the same directory named "dir"
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
#########################################################################################
import argparse
import netCDF4 as nc
import numpy as np

def compare_ncfiles(rootd, ifname):
  ifnames = ifname.split(',')
  if not rootd == '':
    fileA = rootd + '/' + ifnames[0]
    fileB = rootd + '/' + ifnames[1]
  else:
    fileA = ifnames[0]
    fileB = ifnames[1]
  print("Compare the following two netcdf files:")
  print("  File One: " + str(fileA))
  print("  File Two: " + str(fileB))
  same_var_list = np.array([])
  diff_var_list = np.array([])
  diff_var_max  = np.array([])
  with nc.Dataset(fileA) as srcA:
    with nc.Dataset(fileB) as srcB:
      if not srcA.__dict__ == srcB.__dict__:
        print("Note: Global attributes are different between two files")
      for nameA, dimensionA in srcA.dimensions.items():
        dimensionB = srcB.dimensions[nameA]
        dlenA = len(dimensionA)
        dlenB = len(dimensionB)
        if not dlenA == dlenB:
          print("Note: Dimensions are different between two files")
      for nameB, dimensionB in srcB.dimensions.items():
        dimensionA = srcA.dimensions[nameB]
        dlenA = len(dimensionA)
        dlenB = len(dimensionB)
        if not dlenA == dlenB:
          print("Note: Dimensions are different between two files")
      if not len(srcA.variables.items()) == len(srcB.variables.items()):
        print("Note: Numbers of variables are different between two files")
      for vnameA, variableA in srcA.variables.items():
        variableB = srcB.variables[vnameA]
        if not variableA.filters() == variableB.filters():
          print("Note: Attributes are different for the variable "+vnameA+" between two files")
        if not variableA.dimensions == variableB.dimensions:
          print("Note: Dimensions are different for the variable "+vnameA+" between two files")
        dataA = srcA.variables[vnameA][...]
        dataB = srcB.variables[vnameA][...]
        diff  = dataA - dataB
        diff_max = np.nanmax(np.abs(diff))
        if not np.array_equal(dataA, dataB):
          diff_var_list = np.append(diff_var_list,vnameA)
          diff_var_max  = np.append(diff_var_max,diff_max)
        else:
          same_var_list = np.append(same_var_list,vnameA)
  if len(diff_var_list) > 0:
    print("The following variables are different between two files, attached")
    print("please find the variable names and the maximum difference:")
    num = 0
    for varname in diff_var_list:
      print("  ","{:03d}".format(num+1),"   "+varname.ljust(30)+"   ",diff_var_max[num])
      num += 1
  if len(same_var_list) > 0:
    print("The following variables are identical between two files:")
    num = 0
    for varname in same_var_list:
      print("  ","{:03d}".format(num+1),"   "+varname.ljust(30))
      num += 1
 
if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-rd', '--rootd',  help="name of the root directory", default="")
  ap.add_argument('-if', '--ifname', help="input files or their prefix", required=True)
  MyArgs = ap.parse_args()
  compare_ncfiles(MyArgs.rootd, MyArgs.ifname)
