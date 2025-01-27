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

def plot_and_save(crs, hexagons, name, legend_kwds, output_folder, figsize=(10,5), legend=True, cmap='viridis_r', 
             missing_kwds={"color": "lightgrey", "label": "Missing values",},
             bbox_inches="tight",
             ):
    """
    Plots into a figure and saves figure into a file.

    Parameters
    ----------
    crs : 
        an object of the cartopy.crs.Orthographic class, which determines
        the coordinate reference systems (crs) of a given central latitude and 
        longitude.
    hexagons :
        hexagon GeoJSON file
    name : string
        name of column to use.
    legend_kwds : 
        keyword arguments - 'label' will overwrite the auto-generated label.
    title : string
        name of axes.
    figsize : tuple
        size of figure. Default is (10,5)
    legend : boolean
        whether to plot a legend or not
    cmap : string
        the name of a colormap.
    missing_kwds : dictionary
        keyword arguments specifying the colour for missing values. If None, 
        missing values are not plotted.
    bbox_inches : string
        bounding box in inches - will only save the given portion of the figure.
    """
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    # Only plots if there is at least one value that is not null
    if hexagons[name].isnull().all()==False:
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

if __name__ == "__main__":
    plant_type = str(snakemake.config['plant_type'])
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    demand_excel_path = str(snakemake.input.demand_parameters)
    demand_parameters = pd.read_excel(demand_excel_path,index_col='Demand center')
    demand_centers = demand_parameters.index
    transport_methods = ["trucking", "pipeline"]

    # plot LCOH for each hexagon
    # update central coordinates for area considered
    hexagon_bounds = hexagons.geometry.bounds    
    min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
    max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()

    central_lon = (min_lon + max_lon)/2
    central_lat = (min_lat + max_lat)/2

    crs = ccrs.Orthographic(central_longitude = central_lon, central_latitude= central_lat)
    generators = dict(snakemake.config['generators_dict'])

    output_folder = str(snakemake.output)
    check_folder_exists(output_folder)

    for demand_center in demand_centers:
        # plot lowest LC in each location
        plot_and_save(crs, hexagons, f'{demand_center} lowest cost', 
                    {'label':'LC [euros/kg]'}, output_folder)

        
        for transport_method in transport_methods:
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} production cost',
                        {'label':'Production LC [euros/kg]'}, output_folder)

            #%% plot transportation costs
            if plant_type == "Hydrogen":
                if transport_method == "trucking":
                    hexagons[f'{demand_center} total {transport_method} costs'] =\
                        hexagons[f'{demand_center} {transport_method} transport and conversion costs']+\
                            hexagons[f'{demand_center} road construction costs']
                    
                    plot_and_save(crs, hexagons, f'{demand_center} total {transport_method} costs', 
                        {'label':f'{transport_method} cost [euros/kg]'}, output_folder)
                elif transport_method == "pipeline":
                    plot_and_save(crs, hexagons, f'{demand_center} {transport_method} transport and conversion costs', 
                                {'label':f'{transport_method} cost [euros/kg]'}, output_folder)
            elif plant_type == "Ammonia":
                if transport_method == "trucking":
                    hexagons[f'{demand_center} total {transport_method} costs'] =\
                        hexagons[f'{demand_center} {transport_method} transport costs']+\
                            hexagons[f'{demand_center} road construction costs']
            
                    plot_and_save(crs, hexagons, f'{demand_center} total {transport_method} costs', 
                                {'label':f'{transport_method} cost [euros/kg]'}, output_folder)
                elif transport_method == "pipeline":
                    # -- the below doesn't work for ammonia pipeline as it's null. Might need to change nulls to zero 
                    plot_and_save(crs, hexagons, f'{demand_center} {transport_method} transport costs', 
                                {'label':f'{transport_method} cost [euros/kg]'}, output_folder)

            # %% plot total costs
            # -- the below doesn't work for ammonia pipeline as it's null. Might need to change nulls to zero 
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} total cost',
                        {'label':'LC [euros/kg]'}, output_folder)

            # Electrolyzer capacity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} electrolyzer capacity',
                        {'label': 'Size (MW)'}, output_folder)
        
            # Electrolyzer costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} electrolyzer costs',
                        {'label': '$'}, output_folder)

            # H2 storage capacity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} H2 storage capacity',
                        {'label': 'Size (MW)'}, output_folder)
        
            # H2 storage costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} H2 storage costs',
                        {'label': '$'}, output_folder)
            
            # Battery capcity
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} battery capacity', 
                        {'label': 'Size (MW)'}, output_folder)
        
            # battery costs
            plot_and_save(crs, hexagons, f'{demand_center} {transport_method} battery costs',
                        {'label': '$'}, output_folder)

            for generator in generators.keys():
                # generator capacity
                generator = generator.lower()
        
                plot_and_save(crs, hexagons, f'{demand_center} {transport_method} {generator} capacity',
                            {'label': 'Capacity (MW)'}, output_folder)
        
                # generator costs
                plot_and_save(crs, hexagons, f'{demand_center} {transport_method} {generator} costs',
                            {'label': '$'}, output_folder)
        
    # %% plot water costs
    plot_and_save(crs, hexagons, 'Ocean water costs',
                {'label':'Water cost [euros/kg H2]'}, output_folder)

    plt.ticklabel_format(style='plain')

    plot_and_save(crs, hexagons, 'Freshwater costs',
                {'label':'Water cost [euros/kg H2]'}, output_folder)

    plt.ticklabel_format(style='plain') 