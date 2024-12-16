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
from utils import check_folder_exists

def plot_and_save(crs, name, legend_kwds, figsize=(10,5), legend=True, cmap='viridis_r', 
             missing_kwds={"color": "lightgrey", "label": "Missing values",},
             bbox_inches="tight",
             ):
    """
    Plots into a figure and saves figure into a file.

    Parameters
    ----------
    crs : 
        ...
    column : string
        name of column to use.
    legend_kwds : 
        ...
    title : string
        name of axes.
    figsize : tuple
        size of figure. Default is (10,5)
    legend : boolean
        ...
    cmap : string
        ...
    missing_kwds : dictionary
        ...
    bbox_inches : string
        ...
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column = name,
        legend = legend,
        cmap = cmap,
        legend_kwds = legend_kwds,
        missing_kwds = missing_kwds,
    )

    ax.set_title(name)
    fig.savefig(output_folder + f"/{name}.png", bbox_inches=bbox_inches)
    plt.close()

plant_type = "Ammonia" # -- config call
# plant_type = "Hydrogen" # -- config call
hexagons = gpd.read_file("results/hex_cost_components_DJ_2022.geojson")
demand_excel_path = 'parameters/demand_parameters.xlsx'
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
generators = ["Solar", "Wind"] # snakemake config

output_folder = 'plots/DJ_2022'
check_folder_exists(output_folder)

for demand_center in demand_centers:
    # plot lowest LC in each location
    plot_and_save(crs, f'{demand_center} lowest cost', 
                  {'label':'LC [euros/kg]'})

    
    for transport_method in transport_methods:
        plot_and_save(crs, f'{demand_center} {transport_method} production cost',
                      {'label':'Production LC [euros/kg]'})

        #%% plot transportation costs
        if transport_method == "trucking":
            if plant_type == "Hydrogen":
                hexagons[f'{demand_center} total {transport_method} costs'] =\
                    hexagons[f'{demand_center} {transport_method} transport and conversion costs']+\
                        hexagons[f'{demand_center} road construction costs']
                
                plot_and_save(crs, f'{demand_center} total {transport_method} costs', 
                      {'label':f'{transport_method} cost [euros/kg]'})
            elif plant_type == "Ammonia":
                hexagons[f'{demand_center} total {transport_method} costs'] =\
                    hexagons[f'{demand_center} {transport_method} transport costs']+\
                        hexagons[f'{demand_center} road construction costs']
        
                plot_and_save(crs, f'{demand_center} total {transport_method} costs', 
                            {'label':f'{transport_method} cost [euros/kg]'})
        elif transport_methods == "pipeline":
            if plant_type == "Hydrogen":
                plot_and_save(crs, f'{demand_center} {transport_method} transport and conversion costs', 
                            {'label':f'{transport_method} cost [euros/kg]'})
            elif plant_type == "Ammonia":
                plot_and_save(crs, f'{demand_center} {transport_method} transport costs', 
                            {'label':f'{transport_method} cost [euros/kg]'})

        # %% plot total costs
        # -- the below doesn't work for ammonia pipeline as it's null. Might need to change nulls to zero 
        plot_and_save(crs, f'{demand_center} {transport_method} total cost',
                      {'label':'LC [euros/kg]'})

        # Electrolyzer capacity
        plot_and_save(crs, f'{demand_center} {transport_method} electrolyzer capacity',
                      {'label': 'Size (MW)'})
    
        # Electrolyzer costs
        plot_and_save(crs, f'{demand_center} {transport_method} electrolyzer costs',
                      {'label': '$'})

        # H2 storage capacity
        plot_and_save(crs, f'{demand_center} {transport_method} H2 storage capacity',
                      {'label': 'Size (MW)'})
    
        # H2 storage costs
        plot_and_save(crs, f'{demand_center} {transport_method} H2 storage costs',
                      {'label': '$'})
        
        # Battery capcity
        plot_and_save(crs, f'{demand_center} {transport_method} battery capacity', 
                      {'label': 'Size (MW)'})
    
        # battery costs
        plot_and_save(crs, f'{demand_center} {transport_method} battery costs',
                      {'label': '$'})

        for generator in generators:
            # generator capacity
            generator = generator.lower()
    
            plot_and_save(crs, f'{demand_center} {transport_method} {generator} capacity',
                          {'label': 'Capacity (MW)'})
    
            # generator costs
            plot_and_save(crs, f'{demand_center} {transport_method} {generator} costs',
                          {'label': '$'})
    
# %% plot water costs
plot_and_save(crs, 'Ocean water costs',
              {'label':'Water cost [euros/kg H2]'})

plt.ticklabel_format(style='plain')

plot_and_save(crs, 'Freshwater costs',
              {'label':'Water cost [euros/kg H2]'})

plt.ticklabel_format(style='plain') 