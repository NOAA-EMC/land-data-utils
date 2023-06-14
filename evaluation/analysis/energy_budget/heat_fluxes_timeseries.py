#!/usr/bin/env python3
###################################################################################################
## Create time series plots of various heat fluxes and their related variables for the ufs land 
## driver output, examples:
##   1. python heat_fluxes_timeseries.py -i file_name -p 65
##      The code will create time series plots with data stored in files file_name*.nc for the point
##      with the land index 65;
##   2. python heat_fluxes_timeseries.py -i file_name.nc -p 65 -s static_file.nc
##      The code will create time series plots with data stored in file_name.nc for the point with 
##      the land index 65 and the static information can be found in the static file static_file.nc;
##   3. python heat_fluxes_timeseries.py -i file_name.nc -p 65 -t 2011070500,2011070523 
##      The code will create time series plots with data stored in file_name.nc at the point with 
##      the land index 65 during the period between 2011-07-05 00:00 and 2011-07-05 23:00.
## Author: Zhichang Guo and Michael Barlage, contact: Zhichang.Guo@noaa.gov
###################################################################################################
import argparse
import glob
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import timedelta
from datetime import datetime
from netCDF4 import Dataset

def gen_figure(ifname, sfname, ofname, title, lpt, trange):
#define y-unit to x-unit ratio
  ratio   = 0.6 #0.5
  height  = 5
  width   = 14
  rows    = 3
  cols    = 5
  ip      = 1
  cp      = 1004.64
  if not lpt == '':
    ip = int(lpt)
  ip -= 1
  sub_titles = ["Bare", "Under-canopy", "Leaf", "Vegetation", "Bare/Vegetation"]
  info = ''
  if not sfname == '':
    geonc = Dataset(sfname, "r")
    lons = geonc.variables["longitude"][:]
    lats = geonc.variables["latitude"][:]
    styp = geonc.variables["soil_category"][:]
    vtyp = geonc.variables["vegetation_category"][:]
    elev = geonc.variables["elevation"][:]
    info = 'Longitude: '+str(int(lons[ip]*100.)/100)+', Latitude: '+str(int(lats[ip]*100.)/100) + ', Elevation: '
    info = info + str(int(elev[ip]*100)/100) + ', Vegetation type: ' + str(vtyp[ip]) + ', Soil type: ' + str(styp[ip])
  if not 'NC' in ifname[-3:].upper():
    ifname += '*'
    ifnames = glob.glob(ifname)
    ifnames = sorted(ifnames, key=last_14chars)
  else:
    ifnames = [ifname]
  tds = 0
  for f in ifnames:
    datanc = Dataset(f, "r")
    time_size = datanc.dimensions['time'].size
    tds += time_size
    datanc.close()
  hours = tds
  days = tds/24
  dateFmt = '%d'
  if days > 12:
    dateFmt = '%d'
  else:
    dateFmt = '%HZ'
  if len(ifnames) > 1:
    print("First file: ",ifnames[0])
    print("Last file:  ",ifnames[len(ifnames)-1])
  else:
    print("Input file:  ",ifnames[0])
  print("iloc: ",ip)
  dataX = np.array([])
  timeX = np.array([])
  dataY = np.array([])
  tid = 0
  for f in ifnames:
    datanc = Dataset(f, "r")
    tds_f = datanc.dimensions['time'].size
    if tid == 0:
      tgb    = np.array([])
      tgv    = np.array([])
      tleaf  = np.array([])
      tcan   = np.array([])
      t2b    = np.array([])
      t2v    = np.array([])

      chb    = np.array([])
      chv    = np.array([])
      chleaf = np.array([])
      chuc   = np.array([])
      chb2   = np.array([])
      chv2   = np.array([])

      qair   = np.array([])
      pair   = np.array([])
      eair   = np.array([])
      tair   = np.array([])
      rhoair = np.array([])

      fveg   = np.array([])
      shb    = np.array([])
      shc    = np.array([])
      shg    = np.array([])
      shv    = np.array([])

      shb2   = np.array([])
      shc2   = np.array([])
      shg2   = np.array([])
      shv2   = np.array([])

    timetmp  = datanc.variables["time"][:]
    units    = datanc.variables['time'].getncattr('units')
    words    = units.split(' ')
    time_ini = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M:%S')
    tranges  = trange.split(',')
    if len(tranges) < 2 or (tranges[0] == '' and tranges[1] == ''):
      time_start = time_ini + timedelta(seconds=timetmp[0])
      time_stop  = time_ini + timedelta(seconds=timetmp[len(timetmp)-1])
    else:
      time_start = datetime.strptime(tranges[0], '%Y%m%d%H')
      time_stop  = datetime.strptime(tranges[1], '%Y%m%d%H')

    myFmt = mdates.DateFormatter(dateFmt)
    tgbtmp   = datanc.variables["temperature_bare_grd"][...]
    tgvtmp   = datanc.variables["temperature_veg_grd"][...]
    tleaftmp = datanc.variables["temperature_leaf"][...]
    tcantmp  = datanc.variables["temperature_canopy_air"][...]
    t2btmp   = datanc.variables["temperature_bare_2m"][...]
    t2vtmp   = datanc.variables["temperature_veg_2m"][...]

    chbtmp   = datanc.variables["ch_bare_ground"][...]
    chvtmp   = datanc.variables["ch_vegetated"][...]
    chleaftmp= datanc.variables["ch_leaf"][...]
    chuctmp  = datanc.variables["ch_below_canopy"][...]
    chb2tmp  = datanc.variables["ch_bare_ground_2m"][...]
    chv2tmp  = datanc.variables["ch_vegetated_2m"][...]

    qairtmp = datanc.variables["specific_humidity_forcing"][...]
    pairtmp = datanc.variables["surface_pressure_forcing"][...]
    tairtmp = datanc.variables["temperature_forcing"][...]

    fvegtmp = datanc.variables["vegetation_fraction"][...]
    shbtmp  = datanc.variables["sensible_heat_grd_bar"][...]
    shctmp  = datanc.variables["sensible_heat_leaf"][...]
    shgtmp  = datanc.variables["sensible_heat_grd_veg"][...]

    for tid_f in range(tds_f):
      time_cur = time_ini + timedelta(seconds=timetmp[tid_f])
      if time_cur >= time_start and time_cur <= time_stop:
        dataX  = np.append(dataX,float(tid+1))
        timeX  = np.append(timeX,time_cur)
        tgb    = np.append(tgb,tgbtmp[tid_f][ip])
        tgv    = np.append(tgv,tgvtmp[tid_f][ip])
        tleaf  = np.append(tleaf,tleaftmp[tid_f][ip])
        tcan   = np.append(tcan,tcantmp[tid_f][ip])
        t2b    = np.append(t2b,t2btmp[tid_f][ip])
        t2v    = np.append(t2v,t2vtmp[tid_f][ip])
  
        chb    = np.append(chb,chbtmp[tid_f][ip])
        chv    = np.append(chv,chvtmp[tid_f][ip])
        chleaf = np.append(chleaf,chleaftmp[tid_f][ip])
        chuc   = np.append(chuc,chuctmp[tid_f][ip])
        chb2   = np.append(chb2,chb2tmp[tid_f][ip])
        chv2   = np.append(chv2,chv2tmp[tid_f][ip])
  
        qair   = np.append(qair,qairtmp[tid_f][ip])
        pair   = np.append(pair,pairtmp[tid_f][ip])
        tair   = np.append(tair,tairtmp[tid_f][ip])
        etmp   = qairtmp[tid_f][ip]*pairtmp[tid_f][ip]/(0.622+0.378*qairtmp[tid_f][ip])
        eair   = np.append(eair,etmp)
        rair   = (pairtmp[tid_f][ip] - 0.378*etmp)/287.04/tairtmp[tid_f][ip]
        rhoair = np.append(rhoair,rair)

        fveg   = np.append(fveg,fvegtmp[tid_f][ip])
        shb    = np.append(shb,shbtmp[tid_f][ip])
        shc    = np.append(shc,shctmp[tid_f][ip])
        shg    = np.append(shg,shgtmp[tid_f][ip])

        shv    = np.append(shv,shgtmp[tid_f][ip]+shctmp[tid_f][ip]/fvegtmp[tid_f][ip])
        shb2   = np.append(shb2,rair*cp*chbtmp[tid_f][ip]*(tgbtmp[tid_f][ip]-tairtmp[tid_f][ip]))
        shc2   = np.append(shc2,fvegtmp[tid_f][ip]*rair*cp*chleaftmp[tid_f][ip]*(tleaftmp[tid_f][ip]-tcantmp[tid_f][ip]))
        shg2   = np.append(shg2,rair*cp*chuctmp[tid_f][ip]*(tgvtmp[tid_f][ip]-tcantmp[tid_f][ip]))
        shv2   = np.append(shv2,rair*cp*chvtmp[tid_f][ip]*(tcantmp[tid_f][ip]-tairtmp[tid_f][ip]))
        tid += 1
    datanc.close() 
# print(timeX,len(timeX),dataX)

  matplotlib.use("TkAgg")
  fig = plt.figure(figsize=(width,height))
  panel = 0
  for row in range(rows):
    for col in range(cols):
      panel += 1
      ax = fig.add_subplot(rows,cols,panel)
      line1 = 0
      line2 = 0
      line3 = 0
      line4 = 0
      if panel == 1:
        line1 = ax.plot(timeX, shb, label='SHB')
        line2 = ax.plot(timeX, shb2, label='SHB2')
      elif panel == 2:
        line1 = ax.plot(timeX, shg, label='SHG')
        line2 = ax.plot(timeX, shg2, label='SHG2')
      elif panel == 3:
        line1 = ax.plot(timeX, shc, label='SHC')
        line2 = ax.plot(timeX, shc2, label='SHC2')
      elif panel == 4:
        line1 = ax.plot(timeX, shv, label='SHV')
        line2 = ax.plot(timeX, shv2, label='SHV2')
      elif panel == 5:
        line1 = ax.plot(timeX, shb, label='SHB')
        line2 = ax.plot(timeX, shv, label='SHV')
      elif panel == 6:
        line1 = ax.plot(timeX, tair, label='TAIR')
        line2 = ax.plot(timeX, tgb,  label='TGB')
      elif panel == 7:
        line1 = ax.plot(timeX, tcan, label='TCAN')
        line2 = ax.plot(timeX, tgv,  label='TGV')
      elif panel == 8:
        line1 = ax.plot(timeX, tcan, label='TCAN')
        line2 = ax.plot(timeX, tleaf,label='TLEAF')
      elif panel == 9:
        line1 = ax.plot(timeX, tair, label='TAIR')
        line2 = ax.plot(timeX, tcan, label='TCAN')
      elif panel == 10:
        line1 = ax.plot(timeX, t2b,  label='T2B')
        line2 = ax.plot(timeX, t2v,  label='T2V')
        line3 = ax.plot(timeX, tgb,  label='TGB')
        line4 = ax.plot(timeX, tcan, label='TCAN')
      elif panel == 11:
        line1 = ax.plot(timeX, chb,  label='CHB')
      elif panel == 12:
        line1 = ax.plot(timeX, chuc, label='CHUC')
      elif panel == 13:
        line1 = ax.plot(timeX, chleaf,label='CHLEAF')
      elif panel == 14:
        line1 = ax.plot(timeX, chv,  label='CHV')
      elif panel == 15:
        line1 = ax.plot(timeX, chb2, label='CHB2')
        line2 = ax.plot(timeX, chv2, label='CHV2')
      else:
        sys.exit("Error: invalid panel number, should not be here")
      lines = line1
      if not line2 == 0:
        lines += line2
      if not (line3 == 0 or line4 == 0):
        lines += line3 + line4
      ax.xaxis.set_major_formatter(myFmt)
      if days > 12:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=10))
      else:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
      labs = [l.get_label() for l in lines]
      ax.legend(lines, labs, prop={'size': 6})
#get x and y limits
      x_left, x_right = ax.get_xlim()
      y_low, y_high = ax.get_ylim()
#set aspect ratio
      ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
      if row == 0:
        ax.set_title(sub_titles[col], size=10)
      if not row == rows-1:
        ax.set_xticks([])
        ax.set(xlabel=None)
#     else:
#       ax.set(xlabel='Time')
  plt.subplots_adjust(wspace=.022)
  plt.subplots_adjust(hspace=.0)
  plt.subplots_adjust(left = 0.025)
  plt.subplots_adjust(right = 0.99)

  if title.upper() == 'AUTO':
    if info == '':
      plt.suptitle('i: '+str(ip+1), y=0.96, fontsize=13)
    else:
      plt.suptitle('i: '+str(ip+1)+' ('+info+')', y=0.98, fontsize=14)
  elif not title.upper() == 'NONE':
    plt.suptitle(title, y=0.97, fontsize=14)
  if ofname.upper() == 'AUTO':
    plt.savefig('heat_fluxes.'+str(ip+1)+'.png', format='png', dpi=500, bbox_inches='tight',transparent=False)
  elif not ofname.upper() == 'NONE':
    plt.savefig(ofname)
  plt.show()
 
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
def last_14chars(x):
    return(x[-14:])
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-i',   '--input',  help="input files",                      required=True)
    ap.add_argument('-p',   '--point',  help="land point index",                 default='1')
    ap.add_argument('-s',   '--static', help="static input file",                default="")
    ap.add_argument('-o',   '--output', help="output file: none, auto, or name", default="none")
    ap.add_argument('-n',   '--title',  help="figure main title",                default="auto")
    ap.add_argument('-t',   '--trange', help="time range for plotting",          default="")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.static, MyArgs.output, MyArgs.title, MyArgs.point, MyArgs.trange)
