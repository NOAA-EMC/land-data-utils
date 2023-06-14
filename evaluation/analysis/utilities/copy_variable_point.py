#!/usr/bin/env python3
###############################################################################################
## Copy time series of multiple variables from one file to another. The time series are 
## extracted at a grid point over a domain (gridded field data) or a collection of grid points
## (vectorized or aggregated points), either with the grid index or the given longitude and 
## latitude information. examples:
##   1. python copy_variable_point.py -sf source_file.nc -sv tair,tobs -tf target_file.nc -pt
##      207
##      The code will find the grid box with its grid point index 207 and copy its time series
##      of the variables "tair" and "tobs" from the source file "source_file.nc" to the target
##      file "target_file.nc". The target file will be created if it does not exist.
##   2. python copy_variable_point.py -sf src_file.nc -gf static_file.nc -lv lon,lat 
##      -lc 116.69W,36.76 -tf tgt_file.nc -sv temperature,friction_velocity -tv tair,ustar
##      The code will find the nearest grid box with longitude 116.69W and latitude 36.76N. 
##      The longitude and latitude information for the source domain are stored in the static
##      file "static_file.nc" with variable names "lon" and "lat". The time series of the 
##      variables "temperature" and "friction_velocity" in the source file "src_file.nc" will
##      be copied to the target file "tgt_file.nc" with the new variable names "tair" and 
##      "ustar".
##   3. python copy_variable_point.py -sf source_file.nc -tf target_file.nc -sv temperature_soil
##      -tv tsoil2 -pt 207,2
##      The code will find the grid box with its grid point index 207 and copy its time series
##      of the variable "tepmerature_soil" at 2nd layer from the source file "source_file.nc" to
##      the target file "target_file.nc" with the new variable name "tsoil2".
##   4. python copy_variable_point.py -sf source_file.nc -tf target_file.nc -sv temperature_soil
##      -tv tsoil3 -pt 207,99,3 -vc field
##      The code will find the grid box with its i-index 207 and j-index 99 and copy its time 
##      series of the variable "tepmerature_soil" at 3rd layer from the source file 
##      "source_file.nc" to the target file "target_file.nc" with the new variable name "tsoil3".
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
###############################################################################################
from datetime import timedelta
from datetime import datetime
from netCDF4 import Dataset
import argparse
import numpy as np
import sys
import os

def addNote(note, comment):
  if note == None:
    note = comment
  else:
    note += ', '+comment
  print(comment)
  return note

def lonS2F(strLon):
  strLon = strLon.upper()
  strLon = strLon.replace('E','')
  strLon = strLon.replace('N','-')
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

def find_nearest_point(xlon, xlat, lons, lats):
  if xlon < 0.0:
    xlon += 360.
  xlona = xlon - 360.
  xlonb = xlon + 360.
  distanceMin = 9999.9
  index = -1
  found = -1
  dist  = -9.9
  for lid in range(len(lons)):
    if lons[lid] < 0.0:
      lons[lid] += 360.
    distLat = (xlat-lats[lid])*(xlat-lats[lid])
    dist  = (xlon-lons[lid])*(xlon-lons[lid]) + distLat
    distA = (xlona-lons[lid])*(xlona-lons[lid]) + distLat
    distB = (xlonb-lons[lid])*(xlonb-lons[lid]) + distLat
    if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
      distanceMin = min(dist,min(distA,distB))
      index = lid
      found = 1
  if found > 0:
    return index
  else:
    print(strLon+" "+strLat+" "+str(dist)+" "+str(len(lons)))
    sys.exit("the nearest point is not found")

def find_nearest_point2D(xlon, xlat, lons, lats):
  if xlon < 0.0:
    xlon += 360.
  xlona = xlon - 360.
  xlonb = xlon + 360.
  distanceMin = 9999.9
  index = -1
  found = -1
  for yid in range(len(lats)):
    distLat = (xlat-lats[yid])*(xlat-lats[yid])
    for xid in range(len(lons)):
      if lons[xid] < 0.0:
        lons[xid] += 360.
      dist  = (xlon-lons[xid])*(xlon-lons[xid]) + distLat
      distA = (xlona-lons[xid])*(xlona-lons[xid]) + distLat
      distB = (xlonb-lons[xid])*(xlonb-lons[xid]) + distLat
      if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
        distanceMin = min(dist,min(distA,distB))
        i_index = xid
        j_index = yid
        found = 1
  if found > 0:
    return i_index, j_index
  else:
    sys.exit("the nearest point is not found")

def find_nearest_mesh2D(xlon, xlat, lons, lats):
  if xlon < 0.0:
    xlon += 360.
  xlona = xlon - 360.
  xlonb = xlon + 360.
  distanceMin = 9999.9
  index = -1
  found = -1
  for yid in range(len(lats)):
    for xid in range(len(lons)):
      distLat = (xlat-lats[yid,xid])*(xlat-lats[yid,xid])
      if lons[yid,xid] < 0.0:
        lons[yid,xid] += 360.
      dist  = (xlon-lons[yid,xid])*(xlon-lons[yid,xid]) + distLat
      distA = (xlona-lons[yid,xid])*(xlona-lons[yid,xid]) + distLat
      distB = (xlonb-lons[yid,xid])*(xlonb-lons[yid,xid]) + distLat
      if dist <= distanceMin or distA <= distanceMin or distB <= distanceMin:
        distanceMin = min(dist,min(distA,distB))
        i_index = xid
        j_index = yid
        found = 1
  if found > 0:
    return i_index, j_index
  else:
    sys.exit("the nearest point is not found")

def find_tstep_beg(src_time, src_units, tgt_time, tgt_units):
  tstep_beg = -1
  tstep_end = -1
  words = src_units.split(' ')
  src_time_units = words[0]
  src_time_begin = words[2]+' '+words[3]
  src_time_ini = convert_timeStr2timeObj(src_time_begin)
  src_time_beg = advance_dateTime(src_time_ini, src_time[0])
  src_time_end = advance_dateTime(src_time_ini, src_time[len(src_time)-1])
  words = tgt_units.split(' ')
  tgt_time_units = words[0]
  tgt_time_begin = words[2]+' '+words[3]
  tgt_time_ini = convert_timeStr2timeObj(tgt_time_begin)
  tgt_time_beg = advance_dateTime(tgt_time_ini, tgt_time[0])
  tgt_time_end = advance_dateTime(tgt_time_ini, tgt_time[len(tgt_time)-1])
  tstep_beg = None
  for tid in range(len(src_time)):
    time_cur = advance_dateTime(src_time_ini, src_time[tid])
#   if tgt_time_beg == time_cur:
    if tgt_time_beg.year == time_cur.year and tgt_time_beg.month == time_cur.month and tgt_time_beg.day == time_cur.day and tgt_time_beg.hour == time_cur.hour:
      tstep_beg = tid
      break
  if tstep_beg == None:
    sys.exit("Error: time index not found")
  return tstep_beg

def convert_timeStr2timeObj(date_time_str):
# formt = '%d/%m/%y %H:%M:%S' # '18/09/19 01:55:19'
  try:
    formt = '%Y-%m-%d %H:%M:%S' # '2016-01-01 00:50:00'
    date_time_obj = datetime.strptime(date_time_str, formt)
  except:
    formt = '%Y-%m-%d %H:%M' # '2016-01-01 00:50:00'
    date_time_obj = datetime.strptime(date_time_str, formt)
  return date_time_obj

def advance_dateTime(date_time_obj, dt):
  return date_time_obj + timedelta(seconds=dt)

def copy_variable(srcFile, tgtFile, geoFile, srcVar, tgtVar, point, location, 
                  llvname, vector, undef, unitc, mesh):
  rmissing = float(undef)
  if not os.path.isfile(srcFile):
    sys.exit("The source file \""+srcFile+"\" is not found")
  src_grp = Dataset(srcFile, "r") 
  if llvname == '':
    lon_vname = 'longitude'
    lat_vname = 'latitude' 
  else:
    llvns = llvname.split(',')
    lon_vname = llvns[0]
    lat_vname = llvns[1]
  try:
    src_lons = src_grp.variables[lon_vname][:]
    src_lats = src_grp.variables[lat_vname][:]
  except:
    if not os.path.isfile(geoFile):
      sys.exit("The static file \""+geoFile+"\" is not found")
    geo_grp = Dataset(geoFile, "r") 
    src_lons = geo_grp.variables[lon_vname][:]
    src_lats = geo_grp.variables[lat_vname][:]
  note = None
  print("--------------------------------")
  if point == '':
    if location == '':
      if not os.path.isfile(tgtFile):
        sys.exit("no information for finding the point")
      else:
        tgt_grp = Dataset(tgtFile, 'r', format='NETCDF4')
        tgt_lons = tgt_grp.variables['lon'][:]
        tgt_lats = tgt_grp.variables['lat'][:]
        xlon = tgt_lons[0]
        xlat = tgt_lats[0]
        tgt_grp.close()
    else:
      locations = location.split(',')
      xlon = lonS2F(locations[0])
      xlat = latS2F(locations[1])
    if xlon < 0.0:
      xlon += 360.
    print("Find the nearest location for longitude "+str(xlon)+" and latitude "+str(xlat))
    if vector.upper() == 'VECTOR':
      iindex = find_nearest_point(xlon, xlat, src_lons, src_lats)
      jindex = iindex
    else:
      if mesh.upper() == 'TRUE':
        iindex, jindex = find_nearest_mesh2D(xlon, xlat, src_lons, src_lats)
      else:
        iindex, jindex = find_nearest_point2D(xlon, xlat, src_lons, src_lats)
  else:
    lpts = point.split(',')
    iindex = int(lpts[0]) - 1
    jindex = None
    kindex = None
    if len(lpts) > 2:
      kindex = int(lpts[2]) - 1
    if len(lpts) > 1:
      jindex = int(lpts[1]) - 1
    if mesh.upper() == 'TRUE':
      xlon = src_lons[jindex,iindex]
      xlat = src_lats[jindex,iindex]
    else:
      xlon = src_lons[iindex]
      if vector.upper() == 'VECTOR':
        xlat = src_lats[iindex]
      else:
        if jindex == None:
          sys.exit("Error: no information for j-index!")
        xlat = src_lats[jindex]
  if vector.upper() == 'VECTOR':
    comment = "Data copied from the grid with index "+str(iindex+1)
  else:
    comment = "Data copied from the grid with index "+str(iindex+1)+", "+str(jindex+1)
  note = addNote(note, comment)
  if vector.upper() == 'VECTOR':
    comment = "Longitude and latitude for the grid point: "+str(src_lons[iindex])+", "+str(src_lats[iindex])
  else:
    if mesh.upper() == 'TRUE':
      comment = "Longitude and latitude for the grid point: "+str(src_lons[jindex,iindex])+", "+str(src_lats[jindex,iindex])
    else:
      comment = "Longitude and latitude for the grid point: "+str(src_lons[iindex])+", "+str(src_lats[jindex])
  note = addNote(note, comment)

  if tgtVar == '':
    tgtVar = srcVar
  newFileFlag = True
  if not os.path.isfile(tgtFile):
    tgt_grp = Dataset(tgtFile, 'w', format='NETCDF4')
    time = tgt_grp.createDimension('time', None)
    lat  = tgt_grp.createDimension('lat', 1)
    lon  = tgt_grp.createDimension('lon', 1)
    lats = tgt_grp.createVariable('lat', 'f4', ('lat',))
    lats.units = "degrees_north"
    lons = tgt_grp.createVariable('lon', 'f4', ('lon',))
    lons.units = "degrees_east"
    lons[:] = xlon
    lats[:] = xlat
    for vname in src_grp.variables.keys():
      srcVariable = src_grp.variables[vname]
      srcVDims = len(srcVariable.dimensions)
      if 'time' in srcVariable.dimensions and srcVDims == 1:
        x = tgt_grp.createVariable(vname, srcVariable.datatype, srcVariable.dimensions)
        tgt_grp[vname].setncatts(src_grp[vname].__dict__)
        tgt_grp[vname][:] = src_grp[vname][:]  
    tds = len(src_grp.dimensions['time'])
  else:
    newFileFlag = False
    tgt_grp = Dataset(tgtFile, 'a', format='NETCDF4')
    tds = len(tgt_grp.dimensions['time'])
  srcVars = srcVar.split(',')
  tgtVars = tgtVar.split(',')
  unitcs  = unitc.split(',')
  if not len(srcVars) == len(tgtVars):
    sys.exit("Error: number of source and target variables does not match!")
  if not newFileFlag:
    src_times = src_grp.variables['time'][:]
    src_units = src_grp.variables['time'].getncattr('units')
    tgt_times = tgt_grp.variables['time'][:]
    tgt_units = tgt_grp.variables['time'].getncattr('units')
    src_dt    = src_times[1]-src_times[0]
    tgt_dt    = tgt_times[1]-tgt_times[0]
    if tgt_units[:5] == 'hours':
      tgt_dt *= 3600
    if not src_dt == tgt_dt:
      print(src_dt, tgt_dt)
      sys.exit("Error: time interval of the source and target files does not match")
    tstep_beg = find_tstep_beg(src_times, src_units, tgt_times, tgt_units)
  for vid in range(len(srcVars)):
    factor = 1.0
    offset = 0.0
    if vid < len(unitcs):
      unitv = unitcs[vid]
    else:
      unitv = ''
    if not unitv == '':
      if "TIMES" in unitv.upper():
        factor = float(unitv.upper().replace("TIMES",""))
      elif "PLUS" in unitv.upper():
        offset = float(unitv.upper().replace("PLUS",""))
    srcVname= srcVars[vid]
    tgtVname= tgtVars[vid]
    if srcVname in src_grp.variables.keys():
      srcVariable = src_grp.variables[srcVname]
      srcVDims = len(srcVariable.dimensions)
      if 'time' in srcVariable.dimensions and not srcVDims == 1:
        try:
          tgt_grp.createVariable(tgtVname, srcVariable.datatype, ('time', 'lat', 'lon',))
          tgt_grp[tgtVname].setncatts(src_grp[srcVname].__dict__)
        except:
          print("Warning: the variable "+tgtVname+" already exists in the target file and will be overwritten!")
        if not unitv == '':
          tgt_grp[tgtVname].correction = unitv
        src_vars = src_grp.variables[srcVname][:]
        if srcVDims > 4:
          sys.exit("Error: cannot deal with variables with dimensions larger than 4")
        if newFileFlag:
          if srcVDims == 4:
            try:
              if kindex == None:
                sys.exit("Error: no information for k-index!")
              tgt_grp[tgtVname][:] = unit_conversion(src_vars[:,kindex,jindex,iindex],offset,factor,rmissing,unitv)
            except:
              sys.exit("Error: something wrong with reading 4-D variable from the source file")
          elif srcVDims == 3:
            try:
              tgt_grp[tgtVname][:] = unit_conversion(src_vars[:,jindex,iindex], offset, factor, rmissing, unitv)
            except:
              if vector.upper() == 'VECTOR':
                flag = "-pt "+str(iindex+1)+",k-index instead of -lc "+location
              else:
                flag = "-pt "+str(iindex+1)+","+str(jindex+1)+",k-index instead of -lc "+location
              sys.exit("Error: the variable "+tgtVname+" may have vertical dimension, try flag "+flag)
          elif srcVDims == 2:
            tgt_grp[tgtVname][:] = unit_conversion(src_vars[:,iindex], offset, factor, rmissing, unitv)
        else:
          tgt_vars = np.zeros(tds)
          tgt_vars.fill(-9999.)
          tgt_grp[tgtVname][:] = tgt_vars
          k = 0
          for tid in range(tstep_beg,len(src_times)):
            if k < tds:
              if srcVDims == 4:
                try:
                  if kindex == None:
                    sys.exit("Error: no information for k-index!")
                  tgt_vars[k] = unit_conversion_value(src_vars[tid,kindex,jindex,iindex],offset,factor,rmissing,unitv)
                except:
                  sys.exit("Error: something wrong with reading 4-D variable from the source file")
              elif srcVDims == 3:
                try:
                  tgt_vars[k] = unit_conversion_value(src_vars[tid,jindex,iindex],offset,factor,rmissing,unitv)
                except:
                  if vector.upper() == 'VECTOR':
                    flag = "-pt "+str(iindex+1)+",k-index instead of -lc "+location
                  else:
                    flag = "-pt "+str(iindex+1)+","+str(jindex+1)+",k-index instead of -lc "+location
                  sys.exit("Error: the variable "+tgtVname+" may have vertical dimension, try flag "+flag)
              elif srcVDims == 2:
                tgt_vars[k] = unit_conversion_value(src_vars[tid,iindex],offset,factor,rmissing,unitv)
            else:
              break
            k += 1
          tgt_grp[tgtVname][:] = tgt_vars
      else:
        print("Warning: the variable "+srcVname+" skipped as it does not have or only have time dimension!")
    else:
      print("Warning: the variable "+srcVname+" does not exist in the source file!")
  tgt_grp.close()
  src_grp.close()
  print("The program ended normally!")

def not_missing(value, rmissing):
  if value == rmissing or value < rmissing or np.isnan(value):
    return False
  else:
    return True

def unit_conversion(var, offset, factor, rmissing, unitc):
  if not unitc == '':
    for tid in range(len(var)):
      if not_missing(var[tid], rmissing):
        var[tid] = var[tid]*factor + offset
  return var

def unit_conversion_value(var, offset, factor, rmissing, unitc):
  if not unitc == '':
    if not_missing(var, rmissing):
      var = var*factor + offset
  return var
if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-sf', '--srcFile',  help="path to the source file",          required=True)
  ap.add_argument('-tf', '--tgtFile',  help="path to the target file",          required=True)
  ap.add_argument('-sv', '--srcVar',   help="source variable name for copying", required=True)
  ap.add_argument('-gf', '--geoFile',  help="path to the static file",          default='')
  ap.add_argument('-tv', '--tgtVar',   help="target variable name",             default='')
  ap.add_argument('-pt', '--point',    help="location index for copying",       default="")
  ap.add_argument('-lc', '--location', help="longitude/latitude location",      default="")
  ap.add_argument('-lv', '--llvname',  help="longitude/latitude vnames",        default="longitude,latitude")
  ap.add_argument('-vc', '--vector',   help="field or vector",                  default='vector')
  ap.add_argument('-ud', '--undef',    help="value for the missing data",       default='-9999.')
  ap.add_argument('-uc', '--unitc',    help="unit conversion",                  default='')
  ap.add_argument('-mg', '--mesh',     help="mesh grid or not",                 default='FALSE')
  MyArgs = ap.parse_args()
  copy_variable(MyArgs.srcFile, MyArgs.tgtFile, MyArgs.geoFile, MyArgs.srcVar, 
                MyArgs.tgtVar, MyArgs.point, MyArgs.location, MyArgs.llvname, 
                MyArgs.vector, MyArgs.undef, MyArgs.unitc, MyArgs.mesh)
