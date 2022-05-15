import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import math 
from collections import Counter
from numpy import array
import pygrib
import cfgrib
import argparse

#-----CARTOPY---------------------
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

grib= pygrib.open('../gfsanl_4_20160712_0000_000.grb2')

for g in grib:
  #print(g.validDate) 
  time=g.validDate
print(time)

#Variáveis para o plot 
u = grib.select(name='U component of wind', typeOfLevel = 'isobaricInhPa', level = 925)[0]

v= grib.select(name='V component of wind', typeOfLevel = 'isobaricInhPa', level = 925)[0]


# Extract the Brightness mslvalues from the NetCDF 
#(lon mais oeste; lat mais sul; lon mais leste e lat mais norte)     
extent = [-85.0,-60.0,-20.0,10.0]

min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]

ua,lats, lons = u.data(extent[1],extent[3],extent[0],extent[2])
va = v.data(extent[1],extent[3],extent[0],extent[2])[0]

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

#---------- Selecionando as lats e lons para cada var---------

# Correct the longitudes from 0 - 360 
lons_corr = lons - 360

# Calculate the components
Y, X = lats, lons_corr
U, V = ua, va

print(lats,lons)

#---------------------------------------
w=4

S2=ax.streamplot(lons, lats, ua, va,density=3)

#S2 = plt.quiver(lons[::w,::w],lats[::w,::w], ua[::w,::w], va[::w,::w],color='black')

plt.savefig('plot.png')





