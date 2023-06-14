#!/usr/bin/python
######################################################################################
## Count files under a root directory, examples:
##   1. python count_files_rootdir.py -d rootDir
##      Find file counts under the root directory "rootDir"
##   2. python count_files_rootdir.py -d rootDir -r true
##      Find file counts recursively under the root directory "rootDir"
##   3. python count_files_rootdir.py -d rootDir -r true -l 3
##      Find file counts recursively under the root directory "rootDir" with 
##      recursive level 3
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
######################################################################################
from pathlib import Path
import argparse
import os
import sys

def print_indented(result, level):
    print('\t' * level + result)

def traverse_dir(flag, dir):
    iflag = 'not'
    l = os.listdir(dir)
    for d in l:
        if os.path.isdir(dir + d):
            if not d.startswith('.') or (d.startswith('.') and not iflag.upper() == 'NOT'):
                fullpath = os.path.join(dir, d)
                cmd = "find " + " \'" + str(fullpath) + "\'" + " -type f | wc -l "
                result = os.popen(cmd).read()
#               print(str(fullpath)+": "+result,end='')
                print(str(fullpath)+": ",end='')
                print('\033[1;31;43m '+result.rstrip("\n")+' \033[0;0m')

def traverse_dir_recur(flag, dir, max_level, level=0):
    iflag = 'not'
    l = os.listdir(dir)
    for d in l:
        if os.path.isdir(dir + d):
            if not d.startswith('.') or (d.startswith('.') and not iflag.upper() == 'NOT'):
                fullpath = os.path.join(dir, d)
                if max_level == 0 or level < max_level:
                    traverse_dir_recur(flag, dir+  d +"/", max_level, level+1)
                    cmd = "find " + " \'" + str(fullpath) + "\'" + " -type f | wc -l "
                    result = os.popen(cmd).read()
#                   print_indented(str(fullpath)+": "+result, level+1)
#                   print(str(fullpath)+": "+result,end='')
                    print(str(fullpath)+": ",end='')
                    print('\033[1;31;43m '+result.rstrip("\n")+' \033[0;0m')

def find_files(rootdir, rflag, max_level):
    home = os.environ['HOME']
    path = Path(rootdir)
    owner = path.owner()
    print("Owner: ",owner)
    rootdir += '/'
    if rflag == 'false':
        traverse_dir(rflag, rootdir)
    else:
        traverse_dir_recur(rflag, rootdir, max_level)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--rootd', help="name of the root directory", required=True)
    ap.add_argument('-r', '--rflag', help="recurively or not", default="false")
    ap.add_argument('-l', '--level', help="level", type=int, default=0)
    MyArgs = ap.parse_args()
    find_files(MyArgs.rootd, MyArgs.rflag, MyArgs.level)
