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

hexagons = gpd.read_file(str(snakemake.input.hexagons))
demand_excel_path = str(snakemake.input.demand_parameters)
demand_parameters = pd.read_excel(demand_excel_path,index_col='Demand center')
demand_centers = demand_parameters.index
transport_methods = ['pipeline', 'trucking']

# plot LCOH for each hexagon
# update central coordinates for area considered
hexagon_bounds = hexagons.geometry.bounds    
min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()

central_lon = (min_lon + max_lon)/2
central_lat = (min_lat + max_lat)/2

crs = ccrs.Orthographic(central_longitude = central_lon, central_latitude= central_lat)

output_folder = str(snakemake.output)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
for demand_center in demand_centers:
    # plot lowest LCOH in each location
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
    fig.savefig(output_folder + f'/{demand_center} LCOH.png', bbox_inches='tight')
    plt.close()
    
    for transport_method in transport_methods:
        fig = plt.figure(figsize=(10,5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column = f'{demand_center} {transport_method} production cost',
            legend = True,
            cmap = 'viridis_r',
            legend_kwds={'label':'Production LCOH [euros/kg]'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} production cost')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} production cost.png', bbox_inches='tight')
        plt.close()

        #%% plot transportation costs
        fig = plt.figure(figsize=(10,5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons[f'{demand_center} total {transport_method} cost'] =\
            hexagons[f'{demand_center} {transport_method} transport and conversion costs']+hexagons[f'{demand_center} road construction costs']
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column = f'{demand_center} total {transport_method} cost',
            legend = True,
            cmap = 'viridis_r',
            legend_kwds={'label':'{transport_method} cost [euros/kg]'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} transport costs')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} transport cost.png', 
                    bbox_inches='tight')
        plt.close()

        # %% plot total costs
    
        fig = plt.figure(figsize=(10,5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column = f'{demand_center} {transport_method} total cost',
            legend = True,
            cmap = 'viridis_r',
            legend_kwds={'label':'LCOH [euros/kg]'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} LCOH')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} LCOH.png', 
                    bbox_inches='tight')
        plt.close()

        # Electrolyzer capacity
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} electrolyzer capacity',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': 'Size (MW)'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} electrolyzer capacity')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} electrolyzer capacity.png',
                    bbox_inches='tight')
        plt.close()
    
        
        # Electrolyzer costs
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} electrolyzer costs',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': '$'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} electrolyzer costs')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} electrolyzer costs.png',
                    bbox_inches='tight')
        plt.close()

        # H2 storage capacity
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} H2 storage capacity',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': 'Size (MW)'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} H2 storage capacity')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} H2 storage capacity.png', 
                    bbox_inches='tight')
        plt.close()
    
        # H2 storage costs
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} H2 storage costs',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': '$'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} H2 storage costs')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} H2 storage costs.png', 
                    bbox_inches='tight')
        plt.close()

        # Battery capcity
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} battery capacity',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': 'Size (MW)'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} battery capacity')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} battery capacity.png',
                    bbox_inches='tight')
        plt.close()
    
        # battery costs
    
        fig = plt.figure(figsize=(10, 5))
    
        ax = plt.axes(projection=crs)
        ax.set_axis_off()
    
        hexagons.to_crs(crs.proj4_init).plot(
            ax=ax,
            column=f'{demand_center} {transport_method} battery costs',
            legend=True,
            cmap='viridis_r',
            legend_kwds={'label': '$'},
            missing_kwds={
                "color": "lightgrey",
                "label": "Missing values",
            },
        )
        ax.set_title(f'{demand_center} {transport_method} battery costs')
        fig.savefig(output_folder + f'/{demand_center} {transport_method} battery costs.png',
                    bbox_inches='tight')
        plt.close()
        
        for generator in snakemake.config['generators']:
            # generator capacity
            generator = generator.lower()
            fig = plt.figure(figsize=(10, 5))
        
            ax = plt.axes(projection=crs)
            ax.set_axis_off()
        
            hexagons.to_crs(crs.proj4_init).plot(
                ax=ax,
                column=f'{demand_center} {transport_method} {generator} capacity',
                legend=True,
                cmap='viridis_r',
                legend_kwds={'label': 'Capacity (MW)'},
                missing_kwds={
                    "color": "lightgrey",
                    "label": "Missing values",
                },
            )
            ax.set_title(f'{demand_center} {transport_method} {generator} capacity')
            fig.savefig(output_folder + f'/{demand_center} {transport_method} {generator} capacity.png', 
                        bbox_inches='tight')
            plt.close()
    
            # generator costs
        
            fig = plt.figure(figsize=(10, 5))
        
            ax = plt.axes(projection=crs)
            ax.set_axis_off()
        
            hexagons.to_crs(crs.proj4_init).plot(
                ax=ax,
                column=f'{demand_center} {transport_method} {generator} costs',
                legend=True,
                cmap='viridis_r',
                legend_kwds={'label': '$'}, #!!! correct label?
                missing_kwds={
                    "color": "lightgrey",
                    "label": "Missing values",
                },
            )
            ax.set_title(f'{demand_center} {transport_method} {generator} costs ')
            fig.savefig(output_folder + f'/{demand_center} {transport_method} {generator} costs.png', 
                        bbox_inches='tight')
            plt.close()
    
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
plt.ticklabel_format(style='plain')
ax.set_title('Ocean water costs')
fig.savefig(output_folder +'/Ocean water costs.png', bbox_inches='tight')
plt.close()

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
plt.ticklabel_format(style='plain') 
ax.set_title('Freshwater costs')
fig.savefig(output_folder +'/Freshwater costs.png', bbox_inches='tight')
plt.close()


