import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import numpy as np
import scipy
import datetime
import math 
import statistics
import netCDF4
from scipy import stats
from collections import Counter
from numpy import array
import pygrib
import cfgrib
import argparse
from glob import glob
import sys
import os
from numpy import linspace
from numpy import meshgrid
from scipy.ndimage.filters import gaussian_filter
#-----METPY---------------------
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.interpolate import cross_section
from mpl_toolkits.axes_grid1 import make_axes_locatable
from metpy.units import units
#-----CARTOPY---------------------
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader


grib= pygrib.open('../gfsanl_4_20160714_0000_000.grb2')

time='2016-07-14 0000'

#Variaveis 

#Pressure reduced to MSL:Pa cc
pres = grib.select(name='Pressure reduced to MSL')[0]

#U component of wind:m s**-1
u = grib.select(name='U component of wind',   typeOfLevel = 'isobaricInhPa', level = 925)[0]

v= grib.select(name='V component of wind', typeOfLevel = 'isobaricInhPa', level = 925)[0]

#Specific humidity:kg kg**-1
q= grib.select(name='Specific humidity')[0]
#2 metre relative humidity:%
#q2= grib.select(name='2 metre relative humidity')[0]
#t= grib.select(name='Temperature', typeOfLevel = 'isobaricInhPa', level = 925)[0]

# Extract the Brightness mslvalues from the NetCDF 
#(lon mais oeste; lat mais sul; lon mais leste e lat mais norte)     
extent = [-85.0,-60.0,-20.0,10.0]

min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]


#---------- Selecionando as lats e lons para cada var---------
msl, lats, lons = pres.data(min_lat,max_lat,min_lon+360,max_lon+360)

u = u.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]
v = v.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]
  
q =q.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]


#------CALCULO DAS VARIAVEIS------------------------
 
q=q*1000
  
p= (msl/100)

#---------------- CARTOPY --------------------------
fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection=ccrs.PlateCarree())
    
   #Background da imagem 
land = ax.add_feature(cfeature.LAND, facecolor='gray',alpha=0.2)
#ocean = ax.add_feature(cfeature.OCEAN, facecolor='gray')

img_extent = [extent[0], extent[2], extent[1], extent[3]]
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

shapefile = list(shpreader.Reader ('shapefile/ne_10m_admin_1_states_provinces.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(),   edgecolor='gray',facecolor='none', linewidth=0.3)

 # Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='gray',linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='gray', linewidth=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False


#specified interval (n° of lines and columns)
x = linspace(min_lon, max_lon,  u.shape[1])
y = linspace(max_lat, min_lat, v.shape[0]) 

# Create the rectangular grid out of these values
x, y = meshgrid(x, y) 

#-----------------------PLOT--------------------------
## Define de contour interval
min = np.min(q)
max = np.max(q)

interval = 4
levels = np.arange(0,30,interval)


#seleciona valores especificos para colorir diferente
ncolors = list(plt.cm.Greens(np.linspace(0,1,len(levels))))
#intervalo desejado
#ncolors[5] = "white"
#ncolors[6] = "white"
norm = colors.BoundaryNorm(levels,len(levels),ncolors)
cmap = colors.ListedColormap(ncolors,"", len(ncolors))


S1 = ax.contourf(lons,lats,q,cmap=cmap,levels=levels,norm=norm,extend='both')
  
#plt.clabel(S1,inline=1,inline_spacing= 3,fontsize=14,fmt='%1.0f',colors='darkred')
 
# Plot the contour labels
cb2=plt.colorbar(S1,ticks=levels,format='%.0f',pad=0.01, fraction=0.04)

cb2.set_label('(g $kg^{-1}$)',fontsize=14)


S2 = ax.contour(x, y, p, alpha = 0.8, linewidths = 1.5, cmap='copper_r',linestyles='solid')
  
plt.clabel(S2,inline=1,inline_spacing=8,fontsize=12,fmt='%1.0f',colors='darkred')



#-----------Titulos--------------------------------------
plt.title('Umidade específica e PNMM (hPa)', loc='center', fontsize=14,y=1.05)
plt.title('GFS 0.25°', loc='left', fontsize=16)

plt.title('{}'.format(str(pd.to_datetime(str(time), format='%Y-%m-%d %H'))  +  ' UTC'), fontsize=13,loc='right')

plt.savefig('./q/''{}'.format(str(pd.to_datetime(str(time), format='%Y-%m-%d %H')) + "_umidade.png"),dpi=300)
  








