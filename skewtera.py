import datetime 
import scipy.ndimage as ndimage
from siphon.ncss import NCSS
from scipy.ndimage.filters import gaussian_filter
import scipy.integrate as integrate
#Numpy
import numpy as np
from numpy import linspace
from numpy import meshgrid

#Matplotlib
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.cm import get_cmap
import matplotlib.colors as colors
import matplotlib.cm as cm

from netCDF4 import Dataset
from netCDF4 import num2date
import pandas as pd
import xarray as xr

from glob import glob
import sys
import os
#Cartopy
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

#Metpy
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.interpolate import cross_section
from mpl_toolkits.axes_grid1 import make_axes_locatable
from metpy.units import units

import warnings
warnings.filterwarnings("ignore")

#Abrir arquivo 
file =xr.open_dataset("850dia15.nc")
#print(file)

#Período a ser avaliado
date = list(pd.date_range('2016-07-15T00:00:00.00','2016-07-15T02:00:00.00', freq='1H').strftime('%Y-%m-%d %H'))


for idx_tempo in range(len(date)):
  # variavel tempo: string
  ds = file.sel(time=date[idx_tempo]).metpy.parse_cf()
#  print(ds)

