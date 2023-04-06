#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:26:19 2023

@author: Claire Halloran, University of Oxford

Water costs for hydrogen production in each hexagons


"""

import geopandas as gpd
import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt


hexagons = gpd.read_file('Resources/hex_lcoh.geojson')
technology_parameters = "Parameters/technology_parameters.xlsx"
country_excel_path = 'Parameters/country_parameters.xlsx'

water_data = pd.read_excel(technology_parameters,
                            sheet_name='Water',
                            index_col='Parameter'
                            ).squeeze("columns")
country_parameters = pd.read_excel(country_excel_path,
                                    index_col='Country')

#%% water cost for each hexagon for each kg hydrogen produced

h2o_costs_dom_water_bodies = np.empty(len(hexagons))
h2o_costs_ocean = np.empty(len(hexagons))
h2o_costs = np.empty(len(hexagons))

electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
electricity_demand_h2o_ocean_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
water_transport_costs = water_data['Water transport cost (euros/100 km/m3)']
water_spec_cost = water_data['Water specific cost (euros/m3)']
water_demand = water_data['Water demand  (L/kg H2)']

for i in range(len(hexagons)):
    h2o_costs_dom_water_bodies[i] =(water_spec_cost 
                                        + (water_transport_costs/100)*min(hexagons['waterbody_dist'][i],
                                                                          hexagons['waterway_dist'][i]) 
                                        + electricity_demand_h2o_treatment*\
                                            country_parameters.loc[hexagons.country[i],'Electricity price (euros/kWh)']
                                        )*water_demand/1000
    h2o_costs_ocean[i] =(water_spec_cost 
                             + (water_transport_costs/100)*hexagons['ocean_dist'][i] 
                             + electricity_demand_h2o_ocean_treatment*\
                                 country_parameters.loc[hexagons.country[i],'Electricity price (euros/kWh)']
                             )*water_demand/1000
    h2o_costs[i] = min(h2o_costs_dom_water_bodies[i],h2o_costs_ocean[i])

hexagons['Ocean water costs'] = h2o_costs_ocean
hexagons['Freshwater costs'] = h2o_costs_dom_water_bodies
hexagons['Lowest water cost'] = h2o_costs

hexagons.to_file('Resources/hex_water.geojson', driver='GeoJSON')


# %% plot water costs


crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)
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