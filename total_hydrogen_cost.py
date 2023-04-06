#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:44:32 2023

@author: Claire Halloran, University of Oxford

Total hydrogen cost

Bring together all previous data to calculate lowest-cost hydrogen
"""

#%% identify lowest-cost strategy: trucking vs. pipeline

import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
hexagons = gpd.read_file('Resources/hex_water.geojson')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  )

demand_centers = demand_parameters.index
for demand_center in demand_centers:
    hexagons[f'{demand_center} trucking total cost'] =\
        hexagons[f'{demand_center} road construction costs']\
            +hexagons[f'{demand_center} trucking transport and conversion costs']\
                +hexagons[f'{demand_center} trucking production cost']\
                    +hexagons['Lowest water cost']
    hexagons[f'{demand_center} pipeline total cost'] =\
            hexagons[f'{demand_center} pipeline transport and conversion costs']\
                +hexagons[f'{demand_center} pipeline production cost']\
                    +hexagons['Lowest water cost']
                    
    for hexagon in hexagons.index:
        hexagons.loc[hexagon,f'{demand_center} lowest cost'] = np.nanmin(
            [hexagons.loc[hexagon,f'{demand_center} trucking total cost'],
             hexagons.loc[hexagon,f'{demand_center} pipeline total cost']
             ])


# %% plot costs

crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)
for demand_center in demand_centers:
    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} trucking total cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'LCOH [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} trucking LCOH')
    
    crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)
    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} pipeline total cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'LCOH [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} pipeline LCOH')

    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} lowest cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'LCOH [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} LCOH')
