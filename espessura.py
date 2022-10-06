#Importando as LIBS
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
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

#-----CARTOPY---------------------
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader



files= glob("../gfsanl_4_*")

li= sorted(files)
#print(li)

ar=[]
ar2=[]
for idx_tempo in range(len(li)):
  grib= pygrib.open(li[idx_tempo])
  #print(grib)

#Time
  for g in grib:
   
    time=g.validDate
#    print(time)


#Variáveis para o plot 
#press em pa
  pres = grib.select(name='Pressure reduced to MSL')[0]

  sup = grib.select(name='Surface pressure')[0]

#Altura geopot em gpm
  ht_i = grib.select(name='Geopotential Height',   typeOfLevel = 'isobaricInhPa', level = 500)[0]

  ht_f = grib.select(name='Geopotential Height',   typeOfLevel = 'isobaricInhPa', level = 1000)[0]

#U component of wind:m s**-1
  u= grib.select(name='U component of wind', typeOfLevel = 'isobaricInhPa', level = 250)[0]

  v= grib.select(name='V component of wind', typeOfLevel = 'isobaricInhPa', level = 250)[0]


#print(pres)


# Extract the Brightness mslvalues from the NetCDF 
#(lon mais oeste; lat mais sul; lon mais leste e lat mais norte)     
  extent = [-85.0,-60.0,-20.0,10.0]

  min_lon = extent[0]; max_lon = extent[2]; min_lat = extent[1]; max_lat = extent[3]


#---------- Selecionando as lats e lons para cada var---------
  msl, lats, lons = pres.data(min_lat,max_lat,min_lon+360,max_lon+360)

  u = u.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]
  v = v.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]
  
  sup =sup.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]

  ht_i=ht_i.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]
  ht_f=ht_f.data(min_lat,max_lat,min_lon+360,max_lon+360)[0]

#------CALCULO DAS VARIAVEIS----------------------------
  wspd= np.sqrt(u**2 + v**2)
  ht=(ht_i - ht_f)/10
  p=(msl/100)

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

#-----------------------PLOT--------------------------
# Define de contour interval
  var_min = np.min(wspd)
  var_max = np.max(wspd)
  interval = 3
#JAN iniciam-se em 39 m/s
  levels = np.arange(40,var_max,interval)
#------------ Espessura da camada------------------------

  S1 = ax.contour(lons,lats,ht,alpha = 0.6,
  colors='brown',linestyles='dashed')
  
  plt.clabel(S1,inline=1,inline_spacing= 3,fontsize=12,fmt='%1.0f',colors='darkred')
 

#------------------SLP------------------------------------
  S2 = ax.contour(lons,lats,p,colors='black', linestyles='solid', alpha = 0.6) 

  plt.clabel(S2,inline=1,inline_spacing=8,fontsize=12,fmt='%1.0f',colors='black')


#------------Vento em shaded------------------------------
  S3 = ax.contourf(lons,lats, wspd, cmap='BuPu' ,levels=levels,extend='max')

# Plot the contour labels
  cb2=plt.colorbar(S3,ticks=np.arange(40,var_max,interval),format='%.0f',pad=0.01, fraction=0.04)

  cb2.set_label('(m $s^{-1}$)',fontsize=14)
#-----------Titulos--------------------------------------
  plt.title('Vento (m $s^{-1}$) em 250 hPa, Espessura entre 1000 hPa e 500 hPa (dam) e PNMM (hPa)', loc='center', fontsize=14,y=1.05)
  plt.title('GFS 0.25°', loc='left', fontsize=13)

  plt.title('{}'.format(str(pd.to_datetime(str(time), format='%Y-%m-%d %H'))  +  ' UTC'), fontsize=13,loc='right')

  plt.savefig('./esp/''{}'.format(str(pd.to_datetime(str(time), format='%Y-%m-%d %H')) + "_esp.png"))


  print('---------------------------------------')
  print('------ Script started------------------')
  print('File name: ',time)









