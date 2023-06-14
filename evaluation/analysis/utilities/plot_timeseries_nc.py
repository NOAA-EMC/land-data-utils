#!/usr/bin/env python3
###############################################################################################
## plot time series for multiple variables stored in netcdf files, examples:
##   1. python plot_timeseries_nc.py -i file_name.nc -v tair,tobs
##      The code draws time series at the first grid of variables "tair" and "tobs" which are
##      stored in the file "file_name.nc".
##   2. python plot_timeseries_nc.py -i file_name.nc -v tair,tobs -o output.png -p 2,3
##      The code draws time series of variables "tair" and "tobs" which are stored in the file 
##      "file_name.nc". The timeseries are values at grid with 2 and 3 i-index and j-index, 
##      respectively. output.png is the output figure file.
##   3. python plot_timeseries_nc.py -i file_name.nc -v tair -t 2011010100,2011013123 -u n999
##      The code draws time series of the variable "tair" which is stored in the file 
##      "file_name.nc". The timeseries are values at the first grid during the period between
##      2011-01-01 00:00 and 2011-01-31 23:00. The undefined values are labeled with -999.
##   4. python plot_timeseries_nc.py -i file_name.nc -v tair,qair -y L,R -f 1,1000
##      The code draws time series at the first grid of variables "tair" and "qair" which are
##      stored in the file "file_name.nc". tair uses the left y-axis and qair uses the right
##      y-axis. "tair" and "qair" will be factored by 1 and 1000, respectively.
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or Guo.Zhichang@gmail.com
###############################################################################################
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import argparse
from datetime import timedelta
from datetime import datetime
import sys

def plot_time_series(dataX, dataT, dataY, vnames, ofname, yoption, xdate, title):
  width    = 12
  height   = 6
  fig = plt.figure(figsize=(width,height))
  ax = fig.add_subplot(1,1,1)
  yoptions = yoption.split(',')
  dualy = "NO"
  no_ry = 0
  for i in range(len(yoptions)):
    if yoptions[i].upper() == 'R':
      dualy = "YES"
  for vid in range(len(dataY)):
    yside = yoptions[min(vid,len(yoptions)-1)]
    if yside.upper() == 'L':
      if xdate.upper() == 'TRUE':
        line = ax.plot(dataT, dataY[vid], linestyle='solid', label=vnames[vid])
      else:
        line = ax.plot(dataX, dataY[vid], linestyle='solid', label=vnames[vid])
    else:
      if no_ry == 0:
        axr = ax.twinx()
        no_ry = no_ry + 1
      if xdate.upper() == 'TRUE':
        line = axr.plot(dataT, dataY[vid], linestyle='dotted', label=vnames[vid])
      else:
        line = axr.plot(dataX, dataY[vid], linestyle='dotted', label=vnames[vid])
    if vid == 0:
      lines = line
    else:
      lines += line
  diffT = dataT[len(dataT)-1] - dataT[0]
  if diffT.total_seconds() < 400000:
    myLocator = mdates.HourLocator(interval=int(diffT.total_seconds()/18000))
  else:
    myLocator = mdates.DayLocator(interval=int(diffT.total_seconds()/400000))
  if dualy == "YES":
    axr.set(xlabel=None)
  if xdate.upper() == 'TRUE':
    ax.xaxis.set_major_locator(myLocator)
  labs = [l.get_label() for l in lines]
  plt.legend(lines, labs)
  if not title == '':
    plt.suptitle(title,y=0.92,fontsize=14)
  if not ofname == '':
    plt.savefig(ofname)
  else:
    plt.show()
  plt.close('all')

def read_data_mask(datanc, vname, var_mask, mask_id, factor, rmissing):
  vars = datanc.variables[vname][:]
  shape = vars.shape
  dims = len(shape)
  tds = shape[0]
  Y = np.array([])
  if dims == 2 or dims == 3:
    for tid in range(tds):
      avg = cal_masked_average(vars[tid], var_mask, mask_id, rmissing)
      Y = np.append(Y,avg*factor) 
  else:
    sys.exit("Error: cannot deal with dimensions other than 2D and 3D")
  return Y

def cal_masked_average(data, var_mask, mask_id, rmissing):
  newData = data.flatten()
  newMask = var_mask.flatten()
  ma_arr = ma.masked_array(newData, newMask=[(np.isnan(a) or a==rmissing) for a in newData])
  ma_arr = ma.masked_array(ma_arr, mask=[~(a==mask_id) for a in newMask])
  return ma_arr.mean()

def read_data(datanc, vname, iindex, jindex, factor, rmissing):
  vars = datanc.variables[vname][:]
  shape = vars.shape
  dims = len(shape)
  tds = shape[0]
  Y = np.array([])
  if dims == 1:
    for tid in range(tds):
      tmp = vars[tid]
      if abs(tmp-rmissing) > 1.e-8 and not str(tmp)=='--' and not np.isnan(tmp):
        value = tmp*factor
      else:
        value = np.nan
      Y = np.append(Y,value)
  elif dims == 2:
    for tid in range(tds):
      tmp = vars[tid,iindex]
      if abs(tmp-rmissing) > 1.e-8 and not str(tmp)=='--' and not np.isnan(tmp):
        value = tmp*factor
      else:
        value = np.nan
      Y = np.append(Y,value)
  elif dims == 3:
    for tid in range(tds):
      tmp = vars[tid,jindex,iindex]
      if abs(tmp-rmissing) > 1.e-8 and not str(tmp)=='--' and not np.isnan(tmp):
        value = tmp*factor
      else:
        value = np.nan
      Y = np.append(Y,value)
  else:
    sys.exit("Error: cannot deal with dimensions other than 1D, 2D and 3D")
  return Y

def read_var(ifname, gfname, vname, lpt, mask, llvn, rmissing, factor):
  datanc    = nc.Dataset(ifname)
  gf_opened = False
  if gfname == '':
    geonc = datanc
  else:
    geonc = nc.Dataset(gfname)
    gf_opened = True
  mask_name = None
  if not mask == '':
    masks = mask.split("=")
    mask_name = masks[0]
    mask_id   = int(masks[1])
  llvns = llvn.split(',')
  if llvns[0] in geonc.variables.keys() and llvns[1] in geonc.variables.keys():
    lons = geonc.variables[llvns[0]][:]
    lats = geonc.variables[llvns[1]][:]
  elif llvns[0] in datanc.variables.keys() and llvns[1] in datanc.variables.keys():
    lons = datanc.variables[llvns[0]][:]
    lats = datanc.variables[llvns[1]][:]
  else:
    lons = geonc.variables['lon'][:]
    lats = geonc.variables['lat'][:]
  if not mask_name == None:
    if mask_name in geonc.variables.keys():
      var_mask = geonc.variables[mask_name][:]
    elif mask_name in datanc.variables.keys():
      var_mask = datanc.variables[mask_name][:]
  vnames     = vname.split(',') 
  factors    = factor.split(',')
  lpts       = lpt.split(',') 
  tds        = len(datanc.dimensions['time'])
  times      = datanc.variables['time'][:]
  units      = datanc.variables['time'].getncattr('units')
  words      = units.split(' ')
  time_units = words[0]
  try:
    time_ini   = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M:%S')
  except:
    time_ini   = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M')
  iindex     = int(lpts[0]) - 1
  jindex     = int(lpts[min(1,len(lpts)-1)]) - 1

  dataX = np.array([])
  dataT = np.array([])
  dataY = np.array([])
  for tid in range(tds):
    if time_units.upper() == 'DAYS':
      time_cur = time_ini + timedelta(seconds=int(times[tid])*86400)
    elif time_units.upper() == 'HOURS':
      time_cur = time_ini + timedelta(seconds=times[tid]*3600)
    else:
      time_cur = time_ini + timedelta(seconds=times[tid])
    dataX = np.append(dataX,tid)
    dataT = np.append(dataT, time_cur)
  if not mask_name == None:
    for vid in range(len(vnames)):
      vname = vnames[vid]
      f = float(factors[min(vid,len(factors)-1)])
      if vname.startswith('-'):
        vname = vname.replace('-', '')
        f = -1.*f
      Y = read_data_mask(datanc, vname, var_mask, mask_id, f, rmissing)
      if vid == 0:
        dataY = np.array([Y])
      else:
        dataY = np.insert(dataY,len(dataY),Y,0)
  else:
    for vid in range(len(vnames)):
      vname = vnames[vid]
      f = float(factors[min(vid,len(factors)-1)])
      if vname.startswith('-'):
        vname = vname.replace('-', '')
        f = -1.*f
      Y = read_data(datanc, vname, iindex, jindex, f, rmissing)
      if vid == 0:
        dataY = np.array([Y])
      else:
        dataY = np.insert(dataY,len(dataY),Y,0)
  datanc.close()
  if gf_opened:
    geonc.close()
  return dataX, dataT, dataY, vnames

def clip_time(listX, listT, listY, time_beg, time_end):
  newX = []
  newT = []
  newY = []
  for vid in range(len(listY)):
    Y = listY[vid]
    varY = np.array([])
    for tid in range(len(listT)):
      if listT[tid] >= time_beg and listT[tid] <= time_end:
        if vid == 0:
          newX = np.append(newX,listX[tid])
          newT = np.append(newT,listT[tid])
        varY = np.append(varY,Y[tid])
    newY.append(varY)
  return newX, newT, newY

def gen_figure(ifname, gfname, ofname, lpt, vname, mask, llvn, trange, rmissing, yoption, xdate, title, factor):
    rmissing  = rmissing.replace("n","-")
    rmissing  = float(rmissing.replace("N","-"))
    dataX, timeX, dataY, vnames = read_var(ifname, gfname, vname, lpt, mask, llvn, rmissing, factor)
    if not trange == '':
      tranges = trange.split('-')
      time_beg = datetime.strptime(tranges[0], '%Y%m%d%H')
      time_end = datetime.strptime(tranges[1], '%Y%m%d%H')
      dataX, timeX, dataY = clip_time(dataX, timeX, dataY, time_beg, time_end)
      if len(timeX) < 1:
        sys.exit("Warning: no records during this period")
    plot_time_series(dataX, timeX, dataY, vnames, ofname, yoption, xdate, title)

if __name__ == "__main__":
  ap = argparse.ArgumentParser()
  ap.add_argument('-i', '--input',  help="input file path and name",          required=True)
  ap.add_argument('-v', '--vname',  help="variable name for plotting",        required=True)
  ap.add_argument('-g', '--geo',    help="geographic file path and name",     default="")
  ap.add_argument('-o', '--output', help="output file path and name",         default="")
  ap.add_argument('-p', '--point',  help="location index for plotting",       default="1,1,1")
  ap.add_argument('-m', '--mask',   help="mask information",                  default="")
  ap.add_argument('-l', '--llvn',   help="lonitude/latitude variable name",   default="longitude,latitude")
  ap.add_argument('-u', '--undef',  help="missing values",                    default="-9999.0")
  ap.add_argument('-f', '--factor', help="factor",                            default="1")
  ap.add_argument('-y', '--yoption',help="y-axis side option: left or right", default="l")
  ap.add_argument('-t', '--trange', help="time range for plotting",           default="")
  ap.add_argument('-d', '--date',   help="x-axis time date or number",        default="TRUE")
  ap.add_argument('-e', '--title',  help="main title",                        default="")
  MyArgs = ap.parse_args()
  gen_figure(MyArgs.input, MyArgs.geo, MyArgs.output, MyArgs.point, MyArgs.vname, 
             MyArgs.mask, MyArgs.llvn, MyArgs.trange, MyArgs.undef, MyArgs.yoption, 
             MyArgs.date, MyArgs.title, MyArgs.factor)
