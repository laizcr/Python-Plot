#-----------------------------------------------------------------------------------------------------------
import os                            
import requests                     
import time as t                    
import glob    
import datetime 

from datetime import datetime, timedelta, date
from dateutil.relativedelta import relativedelta
    
   
import matplotlib.pyplot as plt          
import cartopy, cartopy.crs as ccrs        
import cartopy.io.shapereader as shpreader 
import numpy as np                         
from netCDF4 import Dataset                         
from matplotlib.patches import Rectangle                                   
import requests 
from itertools import permutations, combinations


import xarray as xr
import pandas as pd
import warnings 
from colour import Color
from datetime import datetime, timedelta, date
from matplotlib.patches import Rectangle 
from matplotlib.colors import LinearSegmentedColormap     
import matplotlib
from netCDF4 import Dataset   
import glob
import requests 
import pandas as pd
from google.cloud import storage, bigquery

warnings.filterwarnings("ignore")
import seaborn as sns

import math
import numpy as np

import sys
sys.path.insert(1, '/mnt/c/scripts/tools/')
from bigquery_bucket_tools import *
from plot_tools import *
from tools import *
#-----------------------------------------------------------------------------------------------------------

date_ref = datetime.strptime(sys.argv[1], '%Y-%m-%d')

date_ = date_ref.date().strftime("%Y-%m-%d")

dirt = '/mnt/c/scripts/relatorios_diarios/ventos_ne/results/'

# Climatologia de vento: modelagem-de-precos.climatologies.era5_wind_monthly_climatology_for_projects

project= 'modelagem-de-precos'
dataset_name = 'climatologies'
table = 'era5_wind_monthly_climatology_for_projects'

bq_client = bigquery.Client(project=project)

list_months = date_range(date_ref, 'day', '-', 20).date().strftime("%m"), date_range(date_ref, 'day', '+', 15).date().strftime("%m")

list_months = [str(int(x)) for x in list_months]

clima_month = query_all_month(bq_client, project, dataset_name, table, list_months)

clima_month['project'] = clima_month['project'].replace('UMR', 'UMA') # parques estão com nomes diferentes

# Dados de vento realizado diário: meteorology-series.era5_data_curated.weighted_wind_daily_data_for_projects
project = 'meteorology-series'
dataset_name = 'era5_data_curated'
table = 'weighted_wind_daily_data_for_projects'

bq_client = bigquery.Client(project=project)

data_ini = date_range(date_ref, 'day', '-', 21).date().strftime("%Y-%m-%d")

data_fim = date_range(date_ref, 'day', '-', 1).date().strftime("%Y-%m-%d")

realizado = query_between_date(bq_client, project, dataset_name, table, data_ini, data_fim, 'date_ref')

realizado['project'] = realizado['project'].replace('UMR', 'UMA')

realizado['month'] = realizado['date_ref'].apply(lambda x: int(pd.to_datetime(x).month))

frame = realizado.merge(clima_month, on=['project', 'month'], how = 'inner', suffixes = ['_realizado', '_clima'])

# Previsão 

# Dados de previsão do ECMWF

# Forecast de vento : forecast-curto-prazo.processed_input_data.ecmwf-processed-data

project= 'forecast-curto-prazo'
dataset_name = 'processed_input_data'
table = 'ecmwf-processed-data'

bq_client = bigquery.Client(project=project)

ec_prev = query_equal_date_list_select(bq_client, project, dataset_name, table, date_, 'time_modelo', ['projeto', 'time_modelo', 'timestamp', 'latitude', 'longitude', 'WS100'])

ec_prev= ec_prev.drop('time_modelo', axis=1).set_index(['timestamp']).groupby(['latitude', 'longitude','projeto',pd.Grouper(freq='1D')]).mean().reset_index().sort_values(by=['timestamp', 'projeto','latitude', 'longitude'], axis=0, ascending=True).rename({'timestamp':'date_ref'},axis=1).reset_index(drop=True)

# Tabela com os pontos de grade do ECMWF mais próximos dos complexos: modelagem-de-precos.metadata.closest_four_points_of_each_complex_era5

project= 'modelagem-de-precos'
dataset_name = 'metadata'
table = 'closets_four_points_of_each_complex_ecmwf_wind_high_resolution'

bq_client = bigquery.Client(project=project)


locs_ec = query_all_date(bq_client, project, dataset_name, table)

locs_ec = locs_ec.rename({'latitude_ecmwf':'latitude','longitude_ecmwf':'longitude'},axis=1)

locs_ec['projeto'] = locs_ec['projeto'].replace('UMR', 'UMA') # parques estão com nomes diferentes

merge_ec = ec_prev.merge(locs_ec ,on=['latitude','longitude','projeto'], how='outer')

merge_ec = merge_ec.drop(merge_ec.columns[[0,1, 5,6]],axis = 1)
# Média ponderada pela distância 

merge_ec = weighted_average(merge_ec, 'WS100', 'dist_km', ['date_ref', 'projeto']).rename(columns={'projeto': 'project'})

merge_ec['month'] = merge_ec['date_ref'].apply(lambda x: int(pd.to_datetime(x).month))

frame_ = merge_ec.merge(clima_month, on=['project', 'month'], how = 'inner', suffixes = ['_prev', '_clima'])

frame['anomaly'] = (frame['ws100_avg_realizado'] / frame['ws100_avg_clima']) * 100

frame_['anomaly'] = (frame_['ws_pond'] / frame_['ws100_avg']) * 100


#Projeto
storage_client = storage.Client(project='modelagem-de-precos')

bucket_name = 'relatorios_meteorologicos' 


# Plot 1 -  UMA FIGURA PARA CADA COMPLEXO 
palette = sns.color_palette("tab10")


num_projetos = len(frame['project'].unique())
col = 3
li = -(-num_projetos // col)
print('Rodando o for ...')


frame['date_ref'] = pd.to_datetime(frame['date_ref'])
frame_['date_ref'] = pd.to_datetime(frame_['date_ref'])

for count, projeto in enumerate(frame['project'].unique()):
    
    fig, ax = plt.subplots(figsize=(12, 8))  
    
    print(projeto)

    ax.set_title(f'Complexo: {projeto}', y=1.13, fontsize=20,weight='bold')

    df_projeto_realizado = frame[frame['project'] == projeto]
    df_projeto_previsao = frame_[frame_['project'] == projeto]

    # Realizado
    s1 = df_projeto_realizado.pivot(index='date_ref', columns='project', values=['anomaly'])
    s1.plot(ax=ax,  linewidth=3, color=palette[count])

    # Previsão
    s2 = df_projeto_previsao.pivot(index='date_ref', columns='project', values=['anomaly'])
    s2.plot(ax=ax, linewidth=3, linestyle='--', color=palette[count])


    ax.set_xlabel('Data',labelpad=-40, fontsize=13,weight='bold')
    ax.set_ylabel('Velocidade do Vento [% Climatologia]', fontsize=13,weight='bold')
    ax.legend([ "Realizado","Previsão"], fontsize=11)

    ax.tick_params(axis='both',which='both', labelsize=16)

    ax.legend(["Realizado","Previsão"],loc='upper center',bbox_to_anchor=(0.50,1.11),ncol=2,fontsize=13) ### 

    plt.yticks(np.arange(50,160, 10), fontsize=15)

    plt.tight_layout()

    # save
    plt.savefig(f'{dirt}{date_}_{projeto}.png', bbox_inches='tight', pad_inches=0.1, dpi=150)
    plt.close() 

    #upp para o bucket 
    upload_blob(storage_client, bucket_name, f'{dirt}{date_}_{projeto}.png', f'relatorios_diarios/ventos_ne/{date_}_{projeto}.png')
