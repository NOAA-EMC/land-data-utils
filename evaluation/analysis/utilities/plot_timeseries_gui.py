#!/usr/bin/python
###############################################################################
## The graphical user interface for displaying time series and their
## stastics (diurnal cycle, daily climatology, and monthly climatology)
##   Useage: python plot_timeseries_gui.py -f filename.nc
## Author: Zhichang Guo, email: Zhichang.Guo@noaa.gov or guo.zhichang@gmail.com
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
import argparse
import netCDF4 as nc
import matplotlib.dates as mdates
import matplotlib
import wx
from datetime import timedelta
from datetime import datetime

class MainFrame(wx.Frame):
  def __init__(self, ifname, point, trange):
    self.width_main  = 870
    self.height_main = 350
    self.ifname      = ifname
    self.point       = point
    self.trange      = trange
    super().__init__(None, title='Main Frame')
#   wx.Frame.__init__(self, None, -1,
#                     'MainFrame', size=(self.width_main, self.height_main))
    self.InitUI()
        
  def InitUI(self):   
    self.cbs = []
           
    pnl = wx.Panel(self)
        
    pnl.SetBackgroundColour('#4f5049')
    vbox = wx.BoxSizer(wx.VERTICAL)

    midPan = wx.Panel(pnl)
    midPan.SetBackgroundColour('#ededed')

    vbox.Add(midPan, wx.ID_ANY, wx.EXPAND | wx.ALL, 20)
    pnl.SetSizer(vbox)

    font0 = wx.Font(14, wx.DECORATIVE, wx.FONTSTYLE_NORMAL, wx.NORMAL)
    font2 = wx.Font(22, wx.DECORATIVE, wx.FONTSTYLE_NORMAL, wx.BOLD)
    self.state  = False
    self.axs    = []
    self.datanc = nc.Dataset(self.ifname)
    self.iindex = None
    self.jindex = None
    self.kindex = None
    if not self.point == '':
      lpts = self.point.split(',')
      self.iindex = int(lpts[0]) - 1
      if len(lpts) >= 2:
        self.jindex = int(lpts[1]) - 1
      if len(lpts) >= 3:
        self.kindex = int(lpts[2]) - 1
    distros = []
    for vname in self.datanc.variables.keys():
      print("Var: ",vname)
      srcVariable = self.datanc.variables[vname]
      srcVDims = len(srcVariable.dimensions)
      if 'time' in srcVariable.dimensions and srcVDims > 1:
        distros = np.append(distros, vname)
    distros = np.append(distros, distros[0])
    distros[0] = 'None'
    self.times = self.datanc.variables['time'][:]
    units      = self.datanc.variables['time'].getncattr('units')
    words      = units.split(' ')
    self.time_units = words[0]
    self.time_ini   = datetime.strptime(words[2]+' '+words[3], '%Y-%m-%d %H:%M:%S')
    if not self.trange == '':
      tranges = self.trange.split('-')
      self.time_beg = datetime.strptime(tranges[0], '%Y%m%d%H')
      self.time_end = datetime.strptime(tranges[1], '%Y%m%d%H')
    self.statist = ['Original Time Series', 'Diurnal Cycle', 'Daily Average',
                    'Daily Climatology', 'Monthly Average', 'Monthly Climatology', 'Default']
    panels  = ['One', 'Two', 'Three', 'Four', 'Five', 'Six', 'All Panels']
    titles  = ['Modeling Variable', 'Observational Variable', 'Additional Variable', 'Time Series Type']
    self.default_value = [ 't2m', 'tobs', 'tair', 'Original Time Series' ]
    self.current_value = [ '',    '',     '',     '']
    offset  = [ 5,                  -8,                       5,                     10 ]
    text_vars    = []
    text_pans    = []
    xpos_main    = 650
    ypos_main    = 5
    self.vds     = 4
    self.pds     = 7
    self.varname = None
    self.stat    = None
    dx           = 180
    dy           = 25
    x0           = 100
    y0           = 86
    text_title = wx.StaticText(midPan, -1, "Select variables for each panel", 
                            (self.width_main/2-200, y0-80), style=wx.ALIGN_CENTER)
    wx.StaticLine(midPan, -1, pos=(x0-90, y0-55), size=(810,20))
    wx.StaticLine(midPan, -1, pos=(x0-90, y0-56), size=(810,20))
    wx.StaticLine(midPan, -1, pos=(x0-90, y0-57), size=(810,20))
    self.current_value = self.default_value
    for i in range(self.pds):
      if i == self.pds-1:
        text = panels[i]
      else:
        text = "Panel "+panels[i]
      text_pan = wx.StaticText(midPan, -1, text+":", (x0-90, y0+i*dy), 
                 style=wx.ALIGN_CENTER)
      for j in range(self.vds):
        self.cell = i*self.vds + j
        if i == 0:
          text_var = wx.StaticText(midPan, -1, titles[j], (x0+j*dx+offset[j]+24, y0-33), 
                     style=wx.ALIGN_CENTER)
        x = x0 + 10 + j*dx
        y = y0 - 10 + i*dy
        if j == self.vds - 1:
          if i == self.pds - 1:
            type_cur = ''
          else:
            type_cur = self.statist[i]
          cb = wx.ComboBox(midPan, choices=self.statist, value=type_cur, pos=(x, y), 
                size=(170,30), style=wx.CB_DROPDOWN|wx.VSCROLL)
#               size=(170,30), style=wx.CB_DROPDOWN | wx.CB_READONLY)
        else:
          if i == self.pds - 1:
            value_cur = ''
          else:
            if self.current_value[j] in distros:
              value_cur = self.current_value[j]
            else:
              value_cur = distros[0]
          cb = wx.ComboBox(midPan, choices=distros, value=value_cur, pos=(x, y), 
                size=(170,30), style=wx.CB_DROPDOWN|wx.VSCROLL)
#               size=(170,30), style=wx.CB_DROPDOWN | wx.CB_READONLY)
        cb.SetFont(wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL, 0, "Courier"))
        self.Bind(wx.EVT_COMBOBOX, self.OnSelect)
        self.cbs = np.append(self.cbs,cb)
        text_vars = np.append(text_vars, text_var)
      text_pans = np.append(text_pans, text_pan)
    for text_pan in text_pans:
      text_pan.SetFont(font0)
    for text_var in text_vars:
      text_var.SetFont(font0)
    text_title.SetFont(font2)
    self.button = wx.Button(midPan, id = wx.ID_ANY, label ="Plot", 
             pos =(self.width_main/2-70, y0+176), size =(100, 42), name ="Plot")
    self.button.SetFont(font2)
    self.button.Bind(wx.EVT_BUTTON, self.onButton)
        
    self.SetSize((self.width_main, self.height_main))
    self.SetPosition(wx.Point(xpos_main, ypos_main))
    self.SetTitle('Plot Time Series GUI')
#   self.Centre()
    self.Show(True)          
        
  def OnSelect(self, e):
    i = e.GetString()
    self.getSelections()

  def getSelections(self):
    for pid in range(self.pds):
      for vid in range(self.vds):
        k = pid*self.vds + vid
        k_last = (self.pds-1)*self.vds + vid
        value_last = self.cbs[k_last].GetValue()
        value_cur  = self.cbs[k].GetValue()
        if value_cur == 'None':
          value_cur = ''
        if not value_last == '':
          if value_last == 'None':
            value_last = ''
          if pid == self.pds-1:
            self.cbs[k].SetValue('')
          else:
            if vid == self.vds-1 and value_last == 'Default':
              self.cbs[k].SetValue(self.statist[pid])
            else:
              self.cbs[k].SetValue(value_last)
          value_cur = value_last
          if pid == self.pds-1:
            self.current_value[vid] = value_cur
        if vid == self.vds-1:
          if pid == 0:
            stat = value_cur
          else:
            stat += ":" + value_cur
        else:
          if vid == 0:
            vname = value_cur
          else:
            vname += ',' + value_cur
      if pid == 0:
        varname = vname
      else:
        varname += ":" + vname
    self.varname = varname
    self.stat = stat

  def onButton(self, event):
    """
    This method is fired when its corresponding button is pressed
    """
    self.getSelections()
    width  = 12
    height = 7
    vnames = self.varname.split(':')
    stats  = self.stat.split(':')
    if self.state:
#     print(plt.get_fignums())
#     plt.close(self.fig)
#     plt.figure().clear()
      self.fig.clf()
#   else:
    self.fig = plt.figure(figsize=(width,height))
    m = 0
    for i in range(3):
      for j in range(2):
        varnames = vnames[m].split(',')
        statp    = stats[m]
        ax     = self.fig.add_subplot(3,2,m+1)
        if not self.state:
          self.axs = np.append(self.axs,ax)
        dataX = np.array([])
        dataT = np.array([])
        tds   = len(self.times)
        tid_beg = None
        tid_end = None
        for tid in range(tds):
          if self.time_units.upper() == 'DAYS':
            time_cur = self.time_ini + timedelta(seconds=self.times[tid]*86400)
          elif self.time_units.upper() == 'HOURS':
            time_cur = self.time_ini + timedelta(seconds=self.times[tid]*3600)
          else:
            time_cur = self.time_ini + timedelta(seconds=self.times[tid])
          if self.trange == '' or (not self.trange == '' and time_cur >= self.time_beg and time_cur <= self.time_end):
            if tid_beg == None:
              tid_beg = tid
            tid_end = tid
            dataX = np.append(dataX, tid)
            dataT = np.append(dataT, time_cur)
        time_dif = dataT[len(dataT)-1] - dataT[0]
        tsecs  = time_dif.total_seconds()
        if tsecs/86400 <= 5:
          myLocator = mdates.HourLocator(interval=int(tsecs/18000))
        else:
          myLocator = mdates.DayLocator(interval=int(tsecs/400000))
        if statp.upper() == 'DIURNAL CYCLE':
          myFmt = mdates.DateFormatter('%HZ')
          myLocator = mdates.HourLocator(interval=3)
        elif statp.upper() == 'DAILY CLIMATOLOGY':
          myFmt = mdates.DateFormatter('%b')
          myLocator = mdates.MonthLocator(bymonthday=15)
        elif statp.upper() == 'MONTHLY CLIMATOLOGY':
          myFmt = mdates.DateFormatter('%b')
          myLocator = mdates.MonthLocator(bymonthday=15)
        else:
          if tsecs/86400 <= 5:
            myFmt = mdates.DateFormatter('%Y-%m-%d')
          else:
            myFmt = mdates.DateFormatter('%m/%Y')
        cnt = 0
        for vid in range(len(varnames)):
          vname = varnames[vid]
          if not vname == '':
            cnt += 1
            vars = self.datanc.variables[vname][:]
            dims = len(vars.shape)
            Y = np.array([])
            if dims == 2:
              for tid in range(tds):
                if tid >= tid_beg and tid <= tid_end:
                  Y = np.append(Y,vars[tid,self.iindex])
            elif dims == 3:
              for tid in range(tds):
                if tid >= tid_beg and tid <= tid_end:
                  Y = np.append(Y,vars[tid,self.jindex,self.iindex])
            if cnt == 1:
              dataY = np.array([Y])
            else:
              dataY = np.insert(dataY,len(dataY),Y,0)
        dataX, dataT, dataY = cal_stat(dataX, dataT, dataY, self.time_ini, statp)
        for vid in range(len(dataY)):
          Y = dataY[vid]
          line = ax.plot(dataT, Y, label=varnames[vid])
          if vid == 0:
            lines = line
          else:
            lines += line
#       for vid in range(len(dataY)):
#         line = ax.plot(dataT, dataY[vid], label=varnames[vid])
        labs = [l.get_label() for l in lines]
        ax.legend(lines, labs)
        ax.xaxis.set_major_locator(myLocator)
        ax.xaxis.set_major_formatter(myFmt)
        m += 1
    self.state = True
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.show()

def cal_stat(dataX, dataT, dataY, time_ini, stat):
  newY = []
  newX = np.array([])
  newT = np.array([])
  if stat.upper() == 'DIURNAL CYCLE':
    hours = 24
    for ihr in range(hours):
      time_cur = time_ini + timedelta(seconds=3600*ihr)
      newX = np.append(newX, ihr)
      newT = np.append(newT, time_cur)
    for vid in range(len(dataY)):
      Y = dataY[vid]
      avg = np.zeros(hours)
      cnt = np.zeros(hours)
      for tid in range(len(dataT)):
        if Y[tid] > -9990 and not np.isnan(Y[tid]):
          avg[dataT[tid].hour] += Y[tid]
          cnt[dataT[tid].hour] += 1.0
      for ihr in range(hours):
        if cnt[ihr] > 0.5:
          avg[ihr] /= cnt[ihr]
        else:
          avg[ihr] = np.nan
      if vid == 0:
        newY = np.array([avg])
      else:
        newY = np.insert(newY,len(newY),avg,0)
  elif stat.upper() == 'DAILY AVERAGE':
    days = (dataT[len(dataT)-1] - dataT[0]).days + 1
    for idy in range(days):
      time_cur = dataT[0] + timedelta(seconds=86400*idy)
      newX = np.append(newX, idy)
      newT = np.append(newT, time_cur)
    for vid in range(len(dataY)):
      Y = dataY[vid]
      avg = np.zeros(days)
      cnt = np.zeros(days)
      for tid in range(len(dataT)):
        if Y[tid] > -9990 and not np.isnan(Y[tid]):
          idy = (dataT[tid]-dataT[0]).days
          avg[idy] += Y[tid]
          cnt[idy] += 1.0
      for idy in range(days):
        if cnt[idy] > 0.5:
          avg[idy] /= cnt[idy]
        else:
          avg[idy] = np.nan
      if vid == 0:
        newY = np.array([avg])
      else:
        newY = np.insert(newY,len(newY),avg,0)
  elif stat.upper() == 'DAILY CLIMATOLOGY':
    year = 1980
    days = 366
    for idy in range(days):
      time_cur = datetime(year=year,month=1,day=1) + timedelta(seconds=86400*idy)
      newX = np.append(newX, idy)
      newT = np.append(newT, time_cur)
    for vid in range(len(dataY)):
      Y = dataY[vid]
      avg = np.zeros(days)
      cnt = np.zeros(days)
      for tid in range(len(dataT)):
        if Y[tid] > -9990 and not np.isnan(Y[tid]):
          month = dataT[tid].month
          day   = dataT[tid].day
          idy = (datetime(year=year,month=month,day=day)-datetime(year=year,month=1,day=1)).days
          avg[idy] += Y[tid]
          cnt[idy] += 1.0
      for idy in range(days):
        if cnt[idy] > 0.5:
          avg[idy] /= cnt[idy]
        else:
          avg[idy] = np.nan
      idy = 59
      if cnt[idy-1] > 0.5 and cnt[idy+1] > 0.5:
        avg[idy] = 0.5*(avg[idy-1]+avg[idy+1])
      else:
        avg[idy] = np.nan
      if vid == 0:
        newY = np.array([avg])
      else:
        newY = np.insert(newY,len(newY),avg,0)
  elif stat.upper() == 'MONTHLY AVERAGE':
    tds = len(dataT)
    months = (dataT[tds-1].year-dataT[0].year)*12 + dataT[tds-1].month - dataT[0].month + 1
    for imo in range(months):
      year  = int((imo+dataT[0].month-1)/12) + dataT[0].year
      month = (imo+dataT[0].month-1)%12 + 1
      time_cur = datetime(year=year,month=month,day=15)
      newX = np.append(newX, imo)
      newT = np.append(newT, time_cur)
    for vid in range(len(dataY)):
      Y = dataY[vid]
      avg = np.zeros(months)
      cnt = np.zeros(months)
      for tid in range(len(dataT)):
        if Y[tid] > -9990 and not np.isnan(Y[tid]):
          imo = (dataT[tid].year-dataT[0].year)*12 + dataT[tid].month - dataT[0].month
          avg[imo] += Y[tid]
          cnt[imo] += 1.0
      for imo in range(months):
        if cnt[imo] > 0.5:
          avg[imo] /= cnt[imo]
        else:
          avg[imo] = np.nan
      if vid == 0:
        newY = np.array([avg])
      else:
        newY = np.insert(newY,len(newY),avg,0)
  elif stat.upper() == 'MONTHLY CLIMATOLOGY':
    tds = len(dataT)
    months = 12
    year  = time_ini.year
    for imo in range(months):
      time_cur = datetime(year=year,month=imo+1,day=15)
      newX = np.append(newX, imo)
      newT = np.append(newT, time_cur)
    for vid in range(len(dataY)):
      Y = dataY[vid]
      avg = np.zeros(months)
      cnt = np.zeros(months)
      for tid in range(len(dataT)):
        if Y[tid] > -9990 and not np.isnan(Y[tid]):
          imo = dataT[tid].month - 1
          avg[imo] += Y[tid]
          cnt[imo] += 1.0
      for imo in range(months):
        if cnt[imo] > 0.5:
          avg[imo] /= cnt[imo]
        else:
          avg[imo] = np.nan
      if vid == 0:
        newY = np.array([avg])
      else:
        newY = np.insert(newY,len(newY),avg,0)
  else:
    return dataX, dataT, dataY
  return newX, newT, newY

def main(ifname, point, trange):
  print("Start app")
  ex = wx.App()
  print("Start main frame")
  MainFrame(ifname, point, trange)
  print("Start main loop")
  ex.MainLoop()    
  print("Done")

if __name__ == '__main__':
  ap = argparse.ArgumentParser()
  ap.add_argument('-if', '--input',   help="input file path and name",                 required=True)
  ap.add_argument('-pt', '--point',   help="location index for plotting",              default="1,1,1")
  ap.add_argument('-tr', '--trange',  help="time range for plotting",                  default="")
  MyArgs = ap.parse_args()
  main(MyArgs.input, MyArgs.point, MyArgs.trange)   
