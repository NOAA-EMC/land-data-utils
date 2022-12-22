#######################################################################################
## Find the nearest grid point index. longitude/latitude are assumed to be functions of
## grid point for the vector data and functions of i/j index respectively for the 2-D 
## field data. Examples:
##   1. python find_nearest_point.py -i file_name.nc -l 60W,20S
##      This code will find the nearest grid point index for the given longitude 60W 
##      and latitude 20S. "longitude" and "latitude" are assumed to be the variable 
##      names in the file "file_name.nc". 
##   2. python find_nearest_point.py -i file_name.nc -l 60,20 -f field -v lon,lat
##      This code will find the i/j index for the grid box which is closest to the
##      location with the given longitude 60E and latitude 20N. "lon" and "lat" are 
##      assumed to be the variable names for the lontitude/latitude in the file 
##      "file_name.nc". 
##   3. python find_nearest_point.py -i static_fields.nc -m vtyp=16 -l 116.69W,36.77
##      This code will find the nearest grid point index with the vegetation type 
##      "vtyp" 16 for the given longitude 116.69W and latitude 36.77N.
##   4. python find_nearest_point.py -i static_fields.nc -m lsmask=1 -l 60,20
##      This code will find the nearest land point with the grid land-sea mask 
##      "lsmask" 1 for the given longitude 60E and and latitude 20N.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
#######################################################################################
import netCDF4 as nc
import argparse
import numpy as np
import os
import sys

def read_var(fname, varname, mask_name, mask_id):
  mask = np.array([])
  varnames = varname.split(',')
  lon_name = varnames[0]
  lat_name = varnames[1]
  ncfile   = nc.Dataset(fname)
  if lon_name in ncfile.variables.keys():
    lons = ncfile.variables[lon_name][:]
  else:
    sys.exit("Error: variable "+lon_name+" cannot be found")
  if lat_name in ncfile.variables.keys():
    lats = ncfile.variables[lat_name][:]
  else:
    sys.exit("Error: variable "+lat_name+" cannot be found")
  if not mask_name == None:
    if mask_name in ncfile.variables.keys():
      mask = ncfile.variables[mask_name][:]
    else:
      print("Warning: mask variable "+mask_name+" cannot be found in the input file, all grids are candidates")
  ncfile.close()
  return lons, lats, mask
    
def lonlatS2F(strLocation):
  if ',' not in strLocation:
    sys.exit("Error: location is not in the format of longitude,latitude")
  else:
    location = strLocation.upper()
    location.replace('E','')
    location.replace('N','')
    locations = location.split(',')
    if 'W' in locations[0]:
      flon = -1.0*float(locations[0].replace('W',''))
    else:
      flon = float(locations[0])
    if 'S' in locations[1]:
      flat = -1.0*float(locations[1].replace('S',''))
    else:
      flat = float(locations[1])
  return flon,flat

def f2s(value):
  return str(int(value*100)/100.)
    
def lonlatF2S(flon,flat):
  flon = int(flon*100)/100.
  if flon > 180:
    flon = 360. - flon
    flon = int(flon*100)/100.
    strLon = str(flon) + 'W'
  elif flon < 0.0:
    strLon = str(-flon) + 'W'
  else:
    strLon = str(flon) + 'E'
  flat = int(flat*100)/100.
  if flat < 0.0:
    strLat = str(-flat) + 'S'
  else:
    strLat = str(flat) + 'N'
  return strLon+","+strLat

def find_nearest_point(strLocation, lons, lats, mask, mask_id):
  xlon, xlat = lonlatS2F(strLocation)
  if xlon < 0.0:
    xlon += 360.
  xlona = xlon - 360.
  xlonb = xlon + 360.
  distanceMin = 9999.9
  lpt_index = -1
  found = -1
  if not len(lons) == len(lats):
    sys.exit("Error: length of longitude and latitude differs")
  maskit = False
  if not mask_id == None and not len(mask) == 0:
    maskit = True
    if not len(lons) == len(mask):
      sys.exit("Error: length of longitude/latitude and mask differs")
  for lid in range(len(lons)):
    if not maskit or (maskit and mask[lid] == mask_id):
      if lons[lid] < 0.0:
        lons[lid] += 360.
      distLat = (xlat-lats[lid])*(xlat-lats[lid])
      dist  = (xlon-lons[lid])*(xlon-lons[lid]) + distLat
      distA = (xlona-lons[lid])*(xlona-lons[lid]) + distLat
      distB = (xlonb-lons[lid])*(xlonb-lons[lid]) + distLat
      if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
        distanceMin = min(dist,min(distA,distB))
        lpt_index = lid
        found = 1
  if found > 0:
    return lpt_index
  else:
    sys.exit("the nearest point is not found")

def find_nearest_point2D(strLocation, lons, lats, mask, mask_id):
  maskit = False
  if not mask_id == None and not len(mask) == 0:
    maskit = True
  xlon, xlat = lonlatS2F(strLocation)
  if xlon < 0.0:
    xlon += 360.
  xlona = xlon - 360.
  xlonb = xlon + 360.
  distanceMin = 9999.9
  x_index = -1
  y_index = -1
  found = -1
  for xid in range(len(lons)):
    for yid in range(len(lats)):
      if not maskit or (maskit and mask[yid,xid] == mask_id):
        if lons[xid] < 0.0:
          lons[xid] += 360.
        distLat = (xlat-lats[yid])*(xlat-lats[yid])
        dist  = (xlon-lons[xid])*(xlon-lons[xid]) + distLat
        distA = (xlona-lons[xid])*(xlona-lons[xid]) + distLat
        distB = (xlonb-lons[xid])*(xlonb-lons[xid]) + distLat
        if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
          distanceMin = min(dist,min(distA,distB))
          x_index = xid
          y_index = yid
          found = 1
  if found > 0:
    return x_index, y_index
  else:
    sys.exit("the nearest point is not found")

def read_and_find(fname, maskit, fov, varname, location):
  mask_name = None
  mask_id   = None
  if not maskit == '':
    words = maskit.split('=')
    mask_name = words[0]
    mask_id = int(words[1])
  lons, lats, mask = read_var(fname, varname, mask_name, mask_id)
  if 'VECTOR' in fov.upper():
    lpt_index = find_nearest_point(location, lons, lats, mask, mask_id)
    print("|-------------------------------------------------------")
    print("| The nearest point index for location ("+location+") is:      "+str(lpt_index+1))
    print("| Its longitude/latitude in the file:  ("+lonlatF2S(lons[lpt_index],lats[lpt_index])+")")
    if not mask_id == None and not len(mask) == 0:
      print("| Its "+mask_name+" in the file is : "+str(mask[lpt_index]))
    print("|-------------------------------------------------------")
  else:
    i_index, j_index = find_nearest_point2D(location, lons, lats, mask, mask_id)
    print("|----------------------------------------------------------")
    print("| The nearest point i/j index for location("+location+"): "+str(i_index+1)+"/"+str(j_index+1))
    print("| Its longitude/latitude in the file: ("+lonlatF2S(lons[i_index],lats[j_index])+")")
    if not mask_id == None and not len(mask) == 0:
      print("| Its "+mask_name+" in the file is : "+str(mask[j_index,i_index]))
    print("|----------------------------------------------------------")

if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-i', '--input',    help="input file path and name", required=True)
  ap.add_argument('-m', '--maskit',   help="mask variable name and the identification number", default='')
  ap.add_argument('-f', '--vector',   help="field or vector", default='vector')
  ap.add_argument('-v', '--variable', help="variable names for longitude and latitude", default="longitude,latitude")
  ap.add_argument('-l', '--location', help="longitude and latitude of the location", required=True)
  MyArgs = ap.parse_args()
  read_and_find(MyArgs.input, MyArgs.maskit, MyArgs.vector, MyArgs.variable, MyArgs.location)
