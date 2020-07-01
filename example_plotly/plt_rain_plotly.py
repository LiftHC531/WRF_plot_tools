#-*- coding: utf-8 -*-
import sys
sys.path.append("..")
from package import wrf_tools, colorbar_tables #local function
from time import time as cpu_time
import os
import copy
import numpy as np
import glob
import Ngl, Nio
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES
import plotly.graph_objs as go
#from mpl_toolkits.basemap import Basemap
import geocat.comp

def get_rain(ff,t0,t1,minus_ini= True):
    rain_exp0 = getvar(ff, "RAINNC", timeidx=t0) #proudced by microphysics Par.
    rain_con0 = getvar(ff, "RAINC", timeidx=t0)  #produced by cumulus Par.
    print("-"*70+"\nVariable: {}, Units: {}\nVariable: {}, Units: {} "\
          .format(rain_exp0.description,rain_exp0.units, \
                  rain_con0.description,rain_con0.units))
    rain_con0 = rain_con0 + rain_exp0
    rain_exp1 = getvar(ff, "RAINNC", timeidx=t1) #proudced by microphysics Par.
    rain_con1 = getvar(ff, "RAINC", timeidx=t1)  #produced by cumulus Par.
    if minus_ini == False: rain_con0 = 0.0
    rain_con0 = (rain_exp1 + rain_con1) - rain_con0
    return rain_con0

def get_title(times):
    time = []
    date = []
    for t in times:
        tt = t[0:19].replace("-",":").replace("T",":").split(":")
        date.append(wrf_tools.Date_string(tt[0],tt[1], \
                    tt[2],tt[3],tt[4],TW_LST=False,plotly=True))
        time.append(tt[0]+tt[1]+tt[2]+"_"+tt[3]+tt[4]); del tt
    LeftString = "Accumulated from {} to {}".format( \
            str(date[0])[:], str(date[1])[:])
    print("\033[44m{} \033[0m".format(LeftString))
    return LeftString



# Get list of of coastline, country and state lon/lat traces
path = "/home/WRF/shapefile/Shimen/"
ff = path+"COUNTY_MOI_1080617.shp" #"TWN_adm2.shp"
f_dam = path+"290---polygon.shp" # All dam
traces_cc = wrf_tools.get_shp_traces(ff) + wrf_tools.get_shp_traces(f_dam,color="blue")
#--------------------------------------------------------------------
start_time = cpu_time()
undef = np.float64(9.96920996839e+36)
#--------------------------------------------------------------------
File = "wrfout_d03_2008-06-15_12:00:00.mean" #wrf_tools.search_wrf_file() 
ff = nc.Dataset(File); #ff = Nio.open_file(File+".nc") #two ways for read
(ntimes, nlev, nlat, nlon, \
times, lat, lon, res_base) = wrf_tools.info(ff)
aim_pt = [nlat//2, nlon//2] #j, i
rg = 5

#tid = wrf_tools.check_plt_time(ntimes, times)
tid = []
tid.append(8); tid.append(56)
rain = get_rain(ff, tid[0], tid[1]); ff.close()
tt = []
for ID in tid:
    tt.append(str(times[ID].values)[0:19])
leftname = get_title(tt)

# [OUTPUT] Grid on destination rectilinear grid (or read the 1D lat and lon from
#          an other .nc file.
newlat1D_rect = np.linspace(lat.min()+0.05, lat.max()-0.05, nlat)
newlon1D_rect = np.linspace(lon.min()+0.05, lon.max()-0.05, nlon)
var_cmax_rect = geocat.comp.rcm2rgrid(lat[:,:], lon[:,:], rain, newlat1D_rect, newlon1D_rect)

""" Plot """
trace = []
trace.append( go.Contour(
       z=var_cmax_rect,
       x=newlon1D_rect,
       y=newlat1D_rect,
       colorscale=colorbar_tables.cwbrain(plotly=True),
       #colorscale= [[0.0, '#313695'], [0.07692307692307693, '#3a67af'], [0.15384615384615385, '#5994c5'], [0.23076923076923078, '#84bbd8'], [0.3076923076923077, '#afdbea'], [0.38461538461538464, '#d8eff5'], [0.46153846153846156, '#d6ffe1'], [0.5384615384615384, '#fef4ac'], [0.6153846153846154, '#fed987'], [0.6923076923076923, '#fdb264'], [0.7692307692307693, '#f78249'], [0.8461538461538461, '#e75435'], [0.9230769230769231, '#cc2727'], [1.0, '#a50026']],
    #zauto=False,  # custom contour levels
    #zmin=0,      # first contour level
    #zmax=50,        # last contour level  => colorscale is centered about 0

       colorbar={
        "borderwidth": 0,
        "outlinewidth": 0,
        "thickness": 15,
        "tickfont": {"size": 14},
        "title": "mm",
                }, #gives your legend some units

       contours={
        "end": 300,
        "showlines": False,
        "size": 20, #this is your contour interval
        "start": 20 
                },

            ))

data = trace + traces_cc
title = leftname

layout = go.Layout(
    title=title,
    showlegend=False,
    hovermode="closest",        # highlight closest point on hover
    xaxis={
        'zeroline':False,
        'showticklabels':False,
        'showline':True,
        'showgrid':False,
        'ticks':'',
        'range':[lon[aim_pt[0],aim_pt[1]]-rg,lon[aim_pt[0],aim_pt[1]]+rg],  # restrict x-axis to range of lon
          },
    yaxis={
        'zeroline':False,
        'showticklabels':False,
        'showline':True,
        'showgrid':False,
        'ticks':'',
        'range':[lat[aim_pt[0],aim_pt[1]]-rg,lat[aim_pt[0],aim_pt[1]]+rg],  # restrict y-axis to range of lat
          },
    autosize=False,
    width=800,
    height=800,
)
fig = go.Figure(data=data, layout=layout)

fig.write_html('accumulated_rain.html', auto_open=False)

end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print(os.path.basename(__file__)+" has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")

