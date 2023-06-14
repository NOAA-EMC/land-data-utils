#!/usr/bin/python
######################################################################################
## Find disk usage under a root directory, examples:
##   1. python du_rootdir.py -d rootDir
##      Find disk usage under the root directory "rootDir"
##   2. python du_rootdir.py -d rootDir -r true
##      Find disk usage recursively under the root directory "rootDir"
##   3. python du_rootdir.py -d rootDir -r true -l 3
##      Find disk usage recursively under the root directory "rootDir" with 
##      recursive level 3
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
######################################################################################
from pathlib import Path
import argparse
import os
import sys

def print_indented(result, level):
    print('\t' * level + result)

def traverse_dir(dir):
    iflag = 'not'
    l = os.listdir(dir)
    for d in l:
        if os.path.isdir(dir + d):
            if not d.startswith('.') or (d.startswith('.') and not iflag.upper() == 'NOT'):
                fullpath = os.path.join(dir, d)
                cmd = "du -sh " + " \'" + str(fullpath) + "\'"
#               print(cmd)
                result = os.popen(cmd).read()
#               print(result)
                print(result,end='')

def traverse_dir_recur(dir, max_level, level=0):
    iflag = 'not'
    l = os.listdir(dir)
    for d in l:
        if os.path.isdir(dir + d):
            if not d.startswith('.') or (d.startswith('.') and not iflag.upper() == 'NOT'):
                fullpath = os.path.join(dir, d)
                if max_level == 0 or level < max_level:
                    traverse_dir_recur(dir + d + "/", max_level, level+1)
                    cmd = "du -sh " + " \'" + str(fullpath) + "\'"
#                   print(cmd)
                    result = os.popen(cmd).read()
#                   print_indented(result.rstrip("\n"), level+1)
                    if level > 0:
                        print('\033[1;31;43m ' + result.rstrip("\n") + ' \033[0;0m')
                    else:
                        print(result.rstrip("\n"))
                    #print(result)

def find_du(rootdir, rflag, max_level):
    home = os.environ['HOME']
    path = Path(rootdir)
    owner = path.owner()
    print("Owner: ",owner)
    rootdir += '/'
    if rflag == 'false':
        traverse_dir(rootdir)
    else:
        traverse_dir_recur(rootdir, max_level)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--rootd', help="name of the root directory", required=True)
    ap.add_argument('-r', '--rflag', help="recurively or not", default="false")
    ap.add_argument('-l', '--level', help="level", type=int, default=0)
    MyArgs = ap.parse_args()
    find_du(MyArgs.rootd, MyArgs.rflag, MyArgs.level)
