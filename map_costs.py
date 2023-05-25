#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:14:48 2023

@author: Claire Halloran, University of Oxford

This script visualizes the spatial cost of hydrogen for each demand center.
"""

import geopandas as gpd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import os

hexagons = gpd.read_file('Resources/hex_total_cost.geojson')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  )

demand_centers = demand_parameters.index

if not os.path.exists('Plots'):
    os.makedirs('Plots')

#%% plot LCOH for each hexagon
# update central coordinates for area considered
crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)
for demand_center in demand_centers:
    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} trucking production cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'Production LCOH [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} trucking production cost')
    plt.savefig(f'Plots/{demand_center} trucking production cost.png',bbox_inches = 'tight', dpi=300)

    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} pipeline production cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'Production LCOH [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} pipeline production cost')
    
    plt.savefig(f'Plots/{demand_center} pipeline production cost.png',bbox_inches = 'tight', dpi=300)


    #%% plot transportation costs
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    
    hexagons[f'{demand_center} total trucking cost'] =\
        hexagons[f'{demand_center} trucking transport and conversion costs']+hexagons[f'{demand_center} road construction costs']
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} total trucking cost',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'Trucking cost [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} trucking transport costs')
    plt.savefig(f'Plots/{demand_center} trucking transport costs.png',bbox_inches = 'tight', dpi=300)

    
    fig = plt.figure(figsize=(10,5))
    
    ax = plt.axes(projection=crs)
    ax.set_axis_off()
    
    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = f'{demand_center} pipeline transport and conversion costs',
        legend = True,
        cmap = 'viridis_r',
        legend_kwds={'label':'Pipeline cost [euros/kg]'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },    
    )
    ax.set_title(f'{demand_center} pipeline transport costs')
    plt.savefig(f'Plots/{demand_center} pipeline transport costs.png',bbox_inches = 'tight', dpi=300)

    # %% plot total costs

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
    plt.savefig(f'Plots/{demand_center} trucking LCOH.png',bbox_inches = 'tight', dpi=300)

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
    plt.savefig(f'Plots/{demand_center} pipeline LCOH.png',bbox_inches = 'tight', dpi=300)

    
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
    plt.savefig(f'Plots/{demand_center} LCOH.png',bbox_inches = 'tight', dpi=300)

# %% plot water costs

fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'Ocean water costs',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Water cost [euros/kg H2]'},
    missing_kwds={
        "color": "lightgrey",
        "label": "Missing values",
    },    
)
ax.set_title('Ocean water costs')
plt.savefig('Plots/ocean water costs.png',bbox_inches = 'tight', dpi=300)

fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'Freshwater costs',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Water cost [euros/kg H2]'},
    missing_kwds={
        "color": "lightgrey",
        "label": "Missing values",
    },    
)
ax.set_title('Freshwater costs')
plt.savefig('Plots/freshwater costs.png',bbox_inches = 'tight', dpi=300)
