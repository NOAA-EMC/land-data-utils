#!/usr/bin/env python3
######################################################################################################
## Copy variables from one netcdf file to another
##     Usage: python copy_variable_nc.py -s source.nc -t target.nc -v varname
##     The vaiable "varname" is copied from the source file "source.nc" to the target file "target.nc"
##     assuming that the dimensions for "varname" exist in the target file.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
######################################################################################################
from netCDF4 import Dataset
import argparse
import numpy as np
import sys
import os

def copy_variable(srcpath, tgtpath, varname):
    if not os.path.isfile(srcpath):
        sys.exit("The source file \""+srcpath+"\" is not found")
    if not os.path.isfile(tgtpath):
        sys.exit("The target file \""+tgtpath+"\" is not found")
    if not tgtpath == srcpath:
        src_rootgrp = Dataset(srcpath, "r")
        tgt_rootgrp = Dataset(tgtpath, "a")
        varnames = varname.split(',')
        for vid in range(len(varnames)):
            if varnames[vid] in src_rootgrp.variables.keys():
                flag = 1
                if varnames[vid] in tgt_rootgrp.variables.keys():
                    print("The variable \""+ varnames[vid] + "\" already exists in the target file \"" + tgtpath + "\", input 'Y' if you want to replace it, otherwise input 'N'")
                    answer = input("Input your choice: \n")
                    if not 'Y' in answer.upper():
                        flag = 0
                        print("Warning: the variable "+varnames[vid]+" is skipped")
                if 1 == flag:
                    for dimobj in src_rootgrp.dimensions:
                        if not dimobj in tgt_rootgrp.dimensions:
                            sys.exit("Error: "+dname+" is not in the dimension list of the target file \""+tgtpath+"\"")
                    srcvar = src_rootgrp.variables[varnames[vid]]
                    if not varnames[vid] in tgt_rootgrp.variables.keys():
                        tgtvar = tgt_rootgrp.createVariable(varnames[vid], srcvar.datatype, srcvar.dimensions)
                        tgt_rootgrp[varnames[vid]].setncatts(src_rootgrp[varnames[vid]].__dict__)
                    tgt_rootgrp[varnames[vid]][:] = src_rootgrp[varnames[vid]][:]
                    print("    The variable \""+varnames[vid]+"\" in \""+srcpath+"\" is copied to the file \""+tgtpath+"\"")
        src_rootgrp.close()
        tgt_rootgrp.close()
    else:
        sys.exit("Warning: the source and target file are same: \""+tgtpath+"\"")
    print("The script ended normally!")
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-s', '--source', help="path to the source file", required=True)
    ap.add_argument('-t', '--target', help="path to the target file", required=True)
    ap.add_argument('-v', '--variable', help="variable name for copying", required=True)
    MyArgs = ap.parse_args()
    copy_variable(MyArgs.source, MyArgs.target, MyArgs.variable)
