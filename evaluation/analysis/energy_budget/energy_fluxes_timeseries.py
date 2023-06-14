#!/usr/bin/env python3
#####################################################################################################
## Create time series plots for various energy balance terms in ufs land driver output, examples:
##   1. python energy_fluxes_timeseries.py.py -i file_name -p 65 -v 4
##      The code will create time series plots with data stored in files file_name*.nc for the point
##      with the land index 65 and the dynamic vegetation option for running the noahmp is 4;
##   2. python energy_fluxes_timeseries.py -i file_name.nc -p 65 -v 4 -s static_file.nc -o output.png
##      The code will create time series plots with data stored in file_name.nc for the point with
##      land index 65, the dynamic vegetation option for running noahmp is 4, the static information 
##      can be found in the static file static_file.nc, and the output figure is output.png;
##   3. python energy_fluxes_timeseries.py -i file_name.nc -p 65 -v 4 -t 2011070500,2011070523
##      The code will create time series plots with data stored in file_name.nc for the point with
##      land index 65, the dynamic vegetation option for running noahmp is 4. The values are plotted
##      for the period from 2011-07-05 00:00 and 2011-07-05 23:00.
## Author: Zhichang Guo and Michael Barlage, contact: Zhichang.Guo@noaa.gov
#####################################################################################################
import argparse
import glob
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import matplotlib.dates as mdates
from datetime import timedelta
from datetime import datetime
from netCDF4 import Dataset

def gen_figure(ifname, sfname, ofname, title, lpt, dveg, trange):
#define y-unit to x-unit ratio
  ratio   = 0.6 #0.5
  height  = 7.4
  width   = 14.6
  rows    = 3
  cols    = 4
  ip      = 1
  if not lpt == '':
    ip = int(lpt)
  ip -= 1
  sub_titles = ["Bare", "Under-canopy", "Leaf", "Vegetation", "Bare", "Under-canopy", "Leaf", "Vegetation", "Grid", "Grid", "Grid", "Grid"]
  info = ''
# print("Read in data")
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
  tranges = trange.split(',')
  dataX = np.array([])
  timeX = np.array([])
  dataY = np.array([])
  tid = 0
  for f in ifnames:
    datanc = Dataset(f, "r")
    tds_f = datanc.dimensions['time'].size
    if tid == 0:
      sag    = np.array([])
      irb    = np.array([])
      shb    = np.array([])
      evb    = np.array([])
      ghb    = np.array([])
      pahb   = np.array([])
      resid1 = np.array([])

      irg    = np.array([])
      shg    = np.array([])
      evg    = np.array([])
      ghv    = np.array([])
      pahg   = np.array([])
      resid2 = np.array([])

      sav    = np.array([])
      irc    = np.array([])
      shc    = np.array([])
      evc    = np.array([])
      tr     = np.array([])
      pahv   = np.array([])
      canhs  = np.array([])
      resid4 = np.array([])

      swveg  = np.array([])
      irveg  = np.array([])
      shveg  = np.array([])
      lhveg  = np.array([])
      ghveg  = np.array([])
      resid5 = np.array([])

      swnet  = np.array([])
      lwnet  = np.array([])
      shf    = np.array([])
      lhf    = np.array([])
      ghf    = np.array([])
      pah    = np.array([])
      resid3 = np.array([])

      swdn   = np.array([])
      swup   = np.array([])
      lwdn   = np.array([])
      lwup   = np.array([])

    timetmp  = datanc.variables["time"][:]
    units    = datanc.variables['time'].getncattr('units')
    words    = units.split(' ')
    time_ini = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M:%S')
    if len(tranges) < 2 or (tranges[0] == '' and tranges[1] == ''):
      time_start = time_ini + timedelta(seconds=timetmp[0])
      time_stop  = time_ini + timedelta(seconds=timetmp[len(timetmp)-1])
    else:
      time_start = datetime.strptime(tranges[0], '%Y%m%d%H')
      time_stop  = datetime.strptime(tranges[1], '%Y%m%d%H')

    myFmt = mdates.DateFormatter(dateFmt)
    if dveg in [ "4", "5", "9", "10" ]:
      fvegtmp = datanc.variables["max_vegetation_frac"][...]
    else:
      fvegtmp = datanc.variables["vegetation_fraction"][...]
    sagtmp  = datanc.variables["sw_absorbed_ground"][...]
    irbtmp  = datanc.variables["lw_absorbed_grd_bare"][...]
    shbtmp  = datanc.variables["sensible_heat_grd_bar"][...]
    evbtmp  = datanc.variables["latent_heat_grd_bare"][...]
    ghbtmp  = datanc.variables["ground_heat_bare"][...]
    pahbtmp = datanc.variables["precip_adv_heat_grd_b"][...]

    irgtmp  = datanc.variables["lw_absorbed_grd_veg"][...]
    shgtmp  = datanc.variables["sensible_heat_grd_veg"][...]
    evgtmp  = datanc.variables["latent_heat_grd_veg"][...]
    ghvtmp  = datanc.variables["ground_heat_veg"][...]
    pahgtmp = datanc.variables["precip_adv_heat_grd_v"][...]

    savtmp   = datanc.variables["sw_absorbed_veg"][...]
    irctmp   = datanc.variables["lw_absorbed_leaf"][...]
    shctmp   = datanc.variables["sensible_heat_leaf"][...]
    evctmp   = datanc.variables["latent_heat_leaf"][...]
    trtmp    = datanc.variables["latent_heat_trans"][...]
    pahvtmp  = datanc.variables["precip_adv_heat_veg"][...]
    canhstmp = datanc.variables["canopy_heat_storage"][...]

    swdntmp  = datanc.variables["downward_shortwave_forcing"][...]
    albtmp   = datanc.variables["albedo_total"][...]
    emistmp  = datanc.variables["emissivity_total"][...]
    lwdntmp  = datanc.variables["downward_longwave_forcing"][...]
    lwuptmp  = datanc.variables["temperature_radiative"][...]
    shftmp   = datanc.variables["sensible_heat_total"][...]
    lhftmp   = datanc.variables["latent_heat_total"][...]
    ghftmp   = datanc.variables["ground_heat_total"][...]
    pahtmp   = datanc.variables["precip_adv_heat_total"][...]

    for tid_f in range(tds_f):
      time_cur = time_ini + timedelta(seconds=timetmp[tid_f])
      if time_cur >= time_start and time_cur <= time_stop:
        dataX  = np.append(dataX,float(tid+1))
        timeX  = np.append(timeX,time_cur)
        fvegt  = fvegtmp[tid_f][ip]
        sag    = np.append(sag,sagtmp[tid_f][ip])
        irb    = np.append(irb,irbtmp[tid_f][ip])
        shb    = np.append(shb,shbtmp[tid_f][ip])
        evb    = np.append(evb,evbtmp[tid_f][ip])
        ghb    = np.append(ghb,ghbtmp[tid_f][ip])
        pahb   = np.append(pahb,pahbtmp[tid_f][ip])
        resid1t= sagtmp[tid_f][ip] - irbtmp[tid_f][ip] - shbtmp[tid_f][ip]
        resid1t= resid1t - evbtmp[tid_f][ip] - ghbtmp[tid_f][ip] + pahbtmp[tid_f][ip]
        resid1 = np.append(resid1,resid1t)

        irg    = np.append(irg,irgtmp[tid_f][ip])
        shg    = np.append(shg,shgtmp[tid_f][ip])
        evg    = np.append(evg,evgtmp[tid_f][ip])
        ghv    = np.append(ghv,ghvtmp[tid_f][ip])
        pahg   = np.append(pahg,pahgtmp[tid_f][ip])
        resid2t= sagtmp[tid_f][ip] - irgtmp[tid_f][ip] - shgtmp[tid_f][ip]
        resid2t= resid2t - evgtmp[tid_f][ip] - ghvtmp[tid_f][ip] + pahgtmp[tid_f][ip]
        resid2 = np.append(resid2,resid2t)

        savt   = savtmp[tid_f][ip]
        irct   = irctmp[tid_f][ip]
        shct   = shctmp[tid_f][ip]
        evct   = evctmp[tid_f][ip]
        trt    = trtmp[tid_f][ip]
        sav    = np.append(sav,savt)
        irc    = np.append(irc,irct)
        shc    = np.append(shc,shct)
        evc    = np.append(evc,evct)
        tr     = np.append(tr,trt)
        pahvt  = pahvtmp[tid_f][ip]
        pahv   = np.append(pahv,pahvt)
        canhs  = np.append(canhs,canhstmp[tid_f][ip])
        resid4t= savt-irct-shct-evct-trt+pahvt-canhstmp[tid_f][ip]
        resid4 = np.append(resid4,resid4t)

        swvegt = savt+sagtmp[tid_f][ip]
        irvegt = irct+irgtmp[tid_f][ip]
        shvegt = shct+shgtmp[tid_f][ip]
        lhvegt = (evct+trt)+evgtmp[tid_f][ip]
        swveg  = np.append(swveg,swvegt)
        irveg  = np.append(irveg,irvegt)
        shveg  = np.append(shveg,shvegt)
        lhveg  = np.append(lhveg,lhvegt)
        ghveg  = np.append(ghveg,ghvtmp[tid_f][ip])
        resid5t= swvegt - irvegt - shvegt - lhvegt - ghvtmp[tid_f][ip]
        resid5t= resid5t + pahgtmp[tid_f][ip] + pahvt-canhstmp[tid_f][ip]
        resid5 = np.append(resid5,resid5t)

        swupt  = albtmp[tid_f][ip]*swdntmp[tid_f][ip]
        swnett = swdntmp[tid_f][ip] - swupt
        swnet  = np.append(swnet,swnett)
        lwdnt  = lwdntmp[tid_f][ip]*emistmp[tid_f][ip]
        lwupt  = 5.67e-8 * emistmp[tid_f][ip] * lwuptmp[tid_f][ip]**4
        ghft   = -ghftmp[tid_f][ip]
        lwnett = fvegtmp[tid_f][ip]*irgtmp[tid_f][ip]
        lwnett+= (1.-fvegtmp[tid_f][ip])*irbtmp[tid_f][ip] + irct
#       lwnett = lwdnt - lwupt
        lwnet  = np.append(lwnet,-lwnett)
        shf    = np.append(shf,shftmp[tid_f][ip])
        lhf    = np.append(lhf,lhftmp[tid_f][ip])
        ghf    = np.append(ghf,ghft)
        pah    = np.append(pah,pahtmp[tid_f][ip])
        resid3t= swnett-lwnett-shftmp[tid_f][ip]-lhftmp[tid_f][ip]-ghft
        resid3t= resid3t - canhstmp[tid_f][ip]+ pahtmp[tid_f][ip]
        resid3 = np.append(resid3,resid3t)

        swdn   = np.append(swdn,swdntmp[tid_f][ip])
        swup   = np.append(swup,swupt)
        lwdn   = np.append(lwdn,lwdnt)
        lwup   = np.append(lwup,lwupt)
        tid += 1
    datanc.close() 

# print("Plot panels")
  matplotlib.use("TkAgg")
  fig = plt.figure(figsize=(width,height))
  panel = 0
  for row in range(rows):
    for col in range(cols):
      panel += 1
#     print("Panel: ", panel)
      ax = fig.add_subplot(rows,cols,panel)
      line = np.array([None, None, None, None, None, None, None, None])
      if panel == 1:
        line[0] = ax.plot(timeX, sag, label='SAG')
        line[1] = ax.plot(timeX, irb, label='IRB')
        line[2] = ax.plot(timeX, shb, label='SHB')
        line[3] = ax.plot(timeX, evb, label='EVB')
        line[4] = ax.plot(timeX, ghb, label='GHB')
        line[5] = ax.plot(timeX, pahb,label='PAHB')
#       line[6] = ax.plot(timeX, resid1, label='RESIDUAL')
      elif panel == 2:
        line[0] = ax.plot(timeX, sag, label='SAG')
        line[1] = ax.plot(timeX, irg, label='IRG')
        line[2] = ax.plot(timeX, shg, label='SHG')
        line[3] = ax.plot(timeX, evg, label='EVG')
        line[4] = ax.plot(timeX, ghv, label='GHV')
        line[5] = ax.plot(timeX, pahg,label='PAHG')
#       line[6] = ax.plot(timeX, resid2, label='RESIDUAL')
      elif panel == 3:
        line[0] = ax.plot(timeX, sav, label='SAV')
        line[1] = ax.plot(timeX, irc, label='IRC')
        line[2] = ax.plot(timeX, shc, label='SHC')
        line[3] = ax.plot(timeX, evc, label='EVC')
        line[4] = ax.plot(timeX, tr,  label='TR')
        line[5] = ax.plot(timeX, pahv,label='PAHV')
        line[6] = ax.plot(timeX, canhs,label='CANHS')
#       line[7] = ax.plot(timeX, resid4, label='RESIDUAL')
      elif panel == 4:
        line[0] = ax.plot(timeX, swveg, label='SWVEG')
        line[1] = ax.plot(timeX, irveg, label='IRVEG')
        line[2] = ax.plot(timeX, shveg, label='SHVEG')
        line[3] = ax.plot(timeX, lhveg, label='LHVEG')
        line[4] = ax.plot(timeX, ghveg, label='GHVEG')
      elif panel == 5:
        line[0] = ax.plot(timeX, resid1,  label='resid1')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0e'))
      elif panel == 6:
        line[0] = ax.plot(timeX, resid2,  label='resid2')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0e'))
      elif panel == 7:
        line[0] = ax.plot(timeX, resid4,  label='resid4')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0e'))
      elif panel == 8:
        line[0] = ax.plot(timeX, resid5,  label='resid5')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0e'))
      elif panel == 9:
        line[0] = ax.plot(timeX, swnet, label='SWNET')
        line[1] = ax.plot(timeX, lwnet, label='LWEM')
        line[2] = ax.plot(timeX, shf,   label='SHF')
        line[3] = ax.plot(timeX, lhf,   label='LHF')
        line[4] = ax.plot(timeX, ghf,   label='GHF')
      elif panel == 10:
        line[0] = ax.plot(timeX, resid3,  label='resid3')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%0.0e'))
      elif panel == 11:
        line[0] = ax.plot(timeX, swdn,  label='SWDN')
        line[1] = ax.plot(timeX, swup , label='SWUP')
        line[2] = ax.plot(timeX, lwdn,  label='LWDN')
        line[3] = ax.plot(timeX, lwup,  label='LWUP')
        line[4] = ax.plot(timeX, shf,   label='SHF')
        line[5] = ax.plot(timeX, lhf,   label='LHF')
        line[6] = ax.plot(timeX, ghf,   label='GHF')
      else:
        ax.axis('off')
#       print("Note: the panel " + str(panel) + " is empty")
      lines = None
      for m in range(len(line)):
        if m == 0:
          lines = line[m]
        else:
          if not None == line[m]:
            lines += line[m]
      if days > 12:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=max(int(days/6),1)))
      else:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
      ax.xaxis.set_major_formatter(myFmt)
      if not lines == None:
        labs = [l.get_label() for l in lines]
        ax.legend(lines, labs, prop={'size': 6})
        ax.set_title(sub_titles[panel-1], size=10)
#get x and y limits
      x_left, x_right = ax.get_xlim()
      y_low, y_high = ax.get_ylim()
#set aspect ratio
      ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
# print("Plot layouts")
  plt.subplots_adjust(wspace=.24)
  plt.subplots_adjust(hspace=.16)
  plt.subplots_adjust(top = 0.95)
  plt.subplots_adjust(bottom = 0.05)

# print("Plot title and output")
  fig.tight_layout(rect=[0, 0.03, 1, 0.95])
  if title.upper() == 'AUTO':
    plt.suptitle('Energy Balances', y=0.97, fontsize=14)
    if info == '':
      plt.suptitle('i: '+str(ip+1), y=0.97, fontsize=13)
    else:
      plt.suptitle('i: '+str(ip+1)+' ('+info+')', y=0.98, fontsize=14)
  elif not title.upper() == 'NONE':
    plt.suptitle(title, y=0.96, fontsize=14)
  if ofname.upper() == 'AUTO':
    plt.savefig('energy_budgets.'+str(ip+1)+'.png', format='png', dpi=500, bbox_inches='tight',transparent=False)
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
    ap.add_argument('-v',   '--dveg',   help="dynamic vegetation option",        required=True)
    ap.add_argument('-s',   '--static', help="static input file",                default="")
    ap.add_argument('-o',   '--output', help="output file: none, auto, or name", default="none")
    ap.add_argument('-n',   '--title',  help="figure main title",                default="auto")
    ap.add_argument('-p',   '--point',  help="land point index",                 default='1')
    ap.add_argument('-t',   '--trange', help="time range for plotting",          default="")
    MyArgs = ap.parse_args()
    gen_figure(MyArgs.input, MyArgs.static, MyArgs.output, MyArgs.title, MyArgs.point, MyArgs.dveg, MyArgs.trange)
