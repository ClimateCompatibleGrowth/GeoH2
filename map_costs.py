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

#hexagons = gpd.read_file('Resources/hex_total_cost.geojson')
hexagons = gpd.read_file('Resources/hex_cost_components.geojson')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,index_col='Demand center')

demand_centers = demand_parameters.index

# plot LCOH for each hexagon
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
    fig.savefig(f'Resources\\{demand_center} trucking production cost.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} pipeline production cost.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} trucking transport cost.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} pipeline transport cost.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} trucking LCOH.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} pipeline LCOH.png', bbox_inches='tight')
    plt.close()

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
    fig.savefig(f'Resources\\{demand_center} LCOH.png', bbox_inches='tight')
    plt.close()

    # Plot production costs by component

    # Solar capacity

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking solar capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking solar capacity')
    fig.savefig(f'Resources\\{demand_center} trucking solar capacity.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline solar capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline solar capacity')
    fig.savefig(f'Resources\\{demand_center} pipeline solar capacity.png', bbox_inches='tight')
    plt.close()

    # Solar costs

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking solar costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking solar costs ')
    fig.savefig(f'Resources\\{demand_center} trucking solar costs.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline solar costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline solar costs')
    fig.savefig(f'Resources\\{demand_center} pipeline solar costs.png', bbox_inches='tight')
    plt.close()

    # Wind capacity 

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking wind capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking wind capacity')
    fig.savefig(f'Resources\\{demand_center} trucking wind capacity.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline wind capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline wind capacity')
    fig.savefig(f'Resources\\{demand_center} pipeline wind capacity.png', bbox_inches='tight')
    plt.close()
    
    # wind costs

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking wind costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking wind costs ')
    fig.savefig(f'Resources\\{demand_center} trucking wind costs.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline wind costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline wind costs')
    fig.savefig(f'Resources\\{demand_center} pipeline wind costs.png', bbox_inches='tight')
    plt.close()

    # Electrolyzer capacity

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking electrolyzer capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking electrolyzer capacity')
    fig.savefig(f'Resources\\{demand_center} trucking electrolyzer capacity.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline electrolyzer capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline electrolyzer capacity')
    fig.savefig(f'Resources\\{demand_center} pipeline electrolyzer capacity.png', bbox_inches='tight')
    plt.close()
    
    # Electrolyzer costs

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking electrolyzer costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking electrolyzer costs')
    fig.savefig(f'Resources\\{demand_center} trucking electrolyzer costs.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline electrolyzer costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline electrolyzer costs')
    fig.savefig(f'Resources\\{demand_center} pipeline electrolyzer costs.png', bbox_inches='tight')
    plt.close()

    # H2 storage capacity

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking H2 storage capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking H2 storage capacity')
    fig.savefig(f'Resources\\{demand_center} trucking H2 storage capacity.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline H2 storage capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline H2 storage capacity')
    fig.savefig(f'Resources\\{demand_center} pipeline H2 storage capacity.png', bbox_inches='tight')
    plt.close()
    
    # H2 storage costs

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking H2 storage costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking H2 storage costs')
    fig.savefig(f'Resources\\{demand_center} trucking H2 storage costs.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline H2 storage costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline H2 storage costs')
    fig.savefig(f'Resources\\{demand_center} pipeline H2 storage costs.png', bbox_inches='tight')
    plt.close()

    # Battery capcity

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking battery storage capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking battery storage capacity')
    fig.savefig(f'Resources\\{demand_center} trucking battery storage capacity.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline battery storage capacity',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': 'Size (MW)'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline battery storage capacity')
    fig.savefig(f'Resources\\{demand_center} pipeline battery storage capacity.png', bbox_inches='tight')
    plt.close()

    # battery costs

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} trucking battery costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} trucking battery costs')
    fig.savefig(f'Resources\\{demand_center} trucking battery costs.png', bbox_inches='tight')
    plt.close()

    fig = plt.figure(figsize=(10, 5))

    ax = plt.axes(projection=crs)
    ax.set_axis_off()

    hexagons.to_crs(crs.proj4_init).plot(
        ax=ax,
        column=f'{demand_center} pipeline battery costs',
        legend=True,
        cmap='viridis_r',
        legend_kwds={'label': '$'},
        missing_kwds={
            "color": "lightgrey",
            "label": "Missing values",
        },
    )
    ax.set_title(f'{demand_center} pipeline battery costs')
    fig.savefig(f'Resources\\{demand_center} pipeline battery costs.png', bbox_inches='tight')
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
fig.savefig(f'Resources\\Ocean water costs.png', bbox_inches='tight')
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
fig.savefig(f'Resources\\Freshwater costs.png', bbox_inches='tight')
plt.close()


