#!/usr/bin/env python3
###############################################################################################
## plot spatial map of variable over stations, examples:
##   1. python plot_station_map.py -if file_name.nc -vn tair
##      The script draws spatial map of the variable "tair" which is stored in the file 
##      "file_name.nc".
##   2. python plot_station_map.py -if file_name.nc -vn tair -of output.png -ft "Fig. 1"
##      The code draws spatial map of variable "tair" which is stored in the file "file_name.nc". 
##      output.png is the output figure file and the figure has the title "Fig. 1".
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
###############################################################################################
import argparse
import os, sys
from os.path import exists
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import numpy as np
from netCDF4 import Dataset

def read_data(ifname, vname, undef, time):
  if not exists(ifname):
    sys.exit("Warning: something is wrong with the input file")
  datanc = Dataset(ifname, "r")
  lons   = datanc.variables['longitude'][:]
  lats   = datanc.variables['latitude'][:]
  var    = datanc.variables[vname][:]
  dims   = len(var.shape)
  datanc.close()
  lpts   = len(lons)
  lon  = []
  lat  = []
  data = []
  rmissing = -9.99E8
  if dims == 1:
    for ip in range(lpts):
      if not lons[ip] == undef and not lats[ip] == undef and not var[ip] == undef:
        lon.append(lons[ip])
        lat.append(lats[ip])
        data.append(var[0][ip])
  elif dims == 2:
    for ip in range(lpts):
      if not lons[ip] == undef and not lats[ip] == undef and not var[time][ip] == undef:
        lon.append(lons[ip])
        lat.append(lats[ip])
        data.append(var[time][ip])
  else:
    sys.exit("Error: cannot deal with variables other than 1-D and 2-D") 
  return lon, lat, data

def lonS2F(strLon):
  strLon = strLon.upper()
  strLon = strLon.replace('E','')
  if 'W' in strLon:
    flon = -1.0*float(strLon.replace('W',''))
  else:
    flon = float(strLon)
  return flon

def latS2F(strLat):
  strLat = strLat.upper()
  strLat = strLat.replace('N','')
  if 'S' in strLat:
    flat = -1.0*float(strLat.replace('S',''))
  else:
    flat = float(strLat)
  return flat

def gen_figure(ifname, vname, ofname, ftitle, fsize, minmax, centered, cmap, domain, undef, time, msize):
  fsize    = fsize.replace('x','X')
  fsizes   = fsize.split('X')
  width    = float(fsizes[0])
  height   = float(fsizes[1])
  domainLL = domain.split(',')
  lonBeg   = lonS2F(domainLL[0])
  lonEnd   = lonS2F(domainLL[1])
  latBeg   = latS2F(domainLL[2])
  latEnd   = latS2F(domainLL[3])
  scaleWE  = 360/(lonEnd-lonBeg)*10
  scaleSN  = 180/(latEnd-latBeg)*10
  scaleMax = 0.2*max(scaleWE,scaleSN)

  # Compute min/max
  lons, lats, var = read_data(ifname, vname, undef, time)
  if not minmax == '':
    minmax = minmax.replace('n','-')
    minmax = minmax.replace('N','-')
    strMinMax = minmax.split(',')
    vmin = float(strMinMax[0])
    vmax = float(strMinMax[1])
  else:
    if centered == "true":
      vmax = np.nanmax(np.abs(var))
      vmin = -vmax
    else:
      vmin = np.nanmin(var)
      vmax = np.nanmax(var)

  matplotlib.use("TkAgg")
  fig    = plt.figure(figsize=(width, height))
# Initialize the plot pointing to the projection
# ------------------------------------------------
  ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
  gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='-')
  gl.top_labels = False
  gl.xlabel_style = {'size': 10, 'color': 'black'}
  gl.ylabel_style = {'size': 10, 'color': 'black'}

  ax.add_feature(cfeature.COASTLINE)

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
  states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')
  ax.add_feature(states_provinces, edgecolor='black')
  ax.add_feature(cfeature.BORDERS)

# Plot scatter data
# ------------------
  sc = ax.scatter(lons, lats, c=var, s=msize*scaleMax, cmap=cmap, 
                  transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

# Plot colorbar
# --------------
  cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.1, fraction=0.06,)

  ax.set_extent([lonBeg, lonEnd, latBeg, latEnd])
  if not ftitle == '':
    ax.set_title(ftitle, pad=10)
  if not ofname == '':
    plt.savefig(ofname, format='png', dpi=500, bbox_inches='tight', transparent=False)
  plt.show()

if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-if', '--ifname',   help="input file path and name",                required=True)
  ap.add_argument('-vn', '--vname',    help="input variable name",                     required=True)
  ap.add_argument('-of', '--ofname',   help="output file path and name",               default="")
  ap.add_argument('-ft', '--ftitle',   help="figure main title",                       default="")
  ap.add_argument('-fs', '--fsize',    help="figure size",                             default="9.5x6.25")
  ap.add_argument('-mm', '--minmax',   help="color bar minimum and maximum limits",    default="")
  ap.add_argument('-dm', '--domain',   help="longitude range for plotting",            default="235.21,293.06,24.74,49.35")
  ap.add_argument('-ud', '--undef',    help="undefined or missing values", type=float, default=-9.99E8)
  ap.add_argument('-ms', '--msize',    help="marker size factor", type=float,          default=1.0)
  ap.add_argument('-ct', '--centered', help="contour color centered or not",           default="true")
  ap.add_argument('-cm', '--cmap',     help="color map",                               default="bwr")
  ap.add_argument('-t',  '--time',     help="time step for plotting", type=int,        default=0)
  MyArgs = ap.parse_args()
  gen_figure(MyArgs.ifname, MyArgs.vname, MyArgs.ofname, MyArgs.ftitle, MyArgs.fsize, MyArgs.minmax, 
             MyArgs.centered, MyArgs.cmap, MyArgs.domain, MyArgs.undef, MyArgs.time, MyArgs.msize)
  print("The script ended normally!")
