#!/usr/bin/env python3
######################################################################################################
## Remove variables in a netcdf file, examples:
##   1. python rm_variables_nc.py -i input.nc -o out.nc -v varname
##      All variables in the input file "input.nc" are copied to the output file "out.nc" except the
##      variable "varname". It is noted that the removal of variables using this script might not 
##      reduce the file size.
##   2. python rm_variables_nc.py -i input.nc -o out.nc -v varnameA,varnameB
##      All variables in the input file "input.nc" are copied to the output file "out.nc" except the
##      variables "varnameA" and "varnameB".
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
######################################################################################################
from netCDF4 import Dataset
import argparse
import sys
import os

def rm_variable(inpath, outpath, varname):
    if not os.path.isfile(inpath):
        sys.exit("The input file \""+inpath+"\" is not found")
    filepath = outpath
    if not outpath == '' or outpath == inpath:
        if os.path.isfile(outpath):
            print("The output file \""+outpath+"\" is existing or same to the input file \""+inpath+"\"")
            answer = input("Input a file name: \n")
            while not answer or answer == inpath or os.path.isfile(answer):
                print("The file \"" + answer + "\" is existing or same to the input file.")
                answer = input("Input a file name: \n")
            filepath = answer
    src_rootgrp = Dataset(inpath, "r")
    tgt_rootgrp = Dataset(filepath, "w")
    varnames = varname.split(',')
    tgt_rootgrp.setncatts(src_rootgrp.__dict__)
    for name, dimension in src_rootgrp.dimensions.items():
        tgt_rootgrp.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    for vname, variable in src_rootgrp.variables.items():
        if vname not in varnames:
            var = tgt_rootgrp.createVariable(vname, variable.datatype, variable.dimensions)
            tgt_rootgrp[vname][:] = src_rootgrp[vname][:]
            tgt_rootgrp[vname].setncatts(src_rootgrp[vname].__dict__)
    src_rootgrp.close()
    tgt_rootgrp.close()
    print("The script ended normally!")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', help="path to the input file", required=True)
    ap.add_argument('-o', '--output', help="path to the output file", required=True)
    ap.add_argument('-v', '--variable', help="variable names", required=True)
    MyArgs = ap.parse_args()
    rm_variable(MyArgs.input, MyArgs.output, MyArgs.variable)
