# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:47:57 2023

@author: Claire Halloran, University of Oxford

Includes code from Nicholas Salmon, University of Oxford, for optimizing
hydrogen plant capacity.

"""

import atlite
import geopandas as gpd
import pypsa
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import p_H2_aux as aux
from functions import CRF
import numpy as np
import logging
import time

logging.basicConfig(level=logging.ERROR)

# !!! add a way to read demand profiles from any locations within a certain distance of hexagon being considered
def optimize_hydrogen_plant(wind_potential, pv_potential, times, demand_profile,
                            wind_max_capacity, pv_max_capacity, 
                            investment_series, basis_fn = None):
    '''
   Optimizes the size of green hydrogen plant components based on renewable potential, hydrogen demand, and investment parameters. 

    Parameters
    ----------
    wind_potential : xarray DataArray
        1D dataarray of per-unit wind potential in hexagon.
    pv_potential : xarray DataArray
        1D dataarray of per-unit solar potential in hexagon.
    times : xarray DataArray
        1D dataarray with timestamps for wind and solar potential.
    country : string
        country in which hexagon is located.
    demand_profile : pandas DataFrame
        hourly dataframe of hydrogen demand in kg.
    investment_series : pandas Series
        interest rate and lifetime information.
    basis_fn : string, optional
        path to basis function for warmstart. The default is None.

    Returns
    -------
    lcoh : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''    
    # ==================================================================================================================
    # Set up network
    # ==================================================================================================================
    # Import a generic network
    n = pypsa.Network(override_component_attrs=aux.create_override_components())

    # Set the time values for the network
    n.set_snapshots(times)

    # Import the design of the H2 plant into the network
    n.import_from_csv_folder("Parameters/Basic_H2_plant")
    
    # Import demand profile
    # Note: All flows are in MW or MWh, conversions for hydrogen done using HHVs. Hydrogen HHV = 39.4 MWh/t
    # hydrogen_demand = pd.read_excel(demand_path,index_col = 0) # Excel file in kg hydrogen, convert to MWh
    n.add('Load',
          'Hydrogen demand',
          bus = 'Hydrogen',
          p_set = demand_profile[0]/1000*39.4,
          )
    
    # ==================================================================================================================
    # Send the weather data to the model
    # ==================================================================================================================
    #!!! need to implement a max capacity based on land availability-- p_nom_max
    n.generators_t.p_max_pu['Wind'] = wind_potential
    n.generators_t.p_max_pu['Solar'] = pv_potential
    
    n.generators.loc['Wind','p_nom_max'] = wind_max_capacity
    n.generators.loc['Solar','p_nom_max'] = pv_max_capacity
    
    
    # ==================================================================================================================
    # Check if the CAPEX input format in Basic_H2_plant is correct, and fix it up if not
    # ==================================================================================================================
    # specify technology-specific and country-specific WACC and lifetime here
    # !!! need to find a capital cost per MWh of hydrogen storage, using what Nick provided for now
    n.generators.loc['Wind','capital_cost'] = n.generators.loc['Wind','capital_cost']\
        * CRF(investment_series['Wind interest rate'], investment_series['Wind lifetime (years)'])
    n.generators.loc['Solar','capital_cost'] = n.generators.loc['Solar','capital_cost']\
        * CRF(investment_series['Solar interest rate'], investment_series['Solar lifetime (years)'])
    for item in [n.links, n.stores]:
        item.capital_cost = item.capital_cost * CRF(investment_series['Plant interest rate'],investment_series['Plant lifetime (years)'])

    # ==================================================================================================================
    # Solve the model
    # ==================================================================================================================
    solver = 'gurobi'
    if basis_fn is None:
        n.lopf(solver_name=solver,
               solver_options = {'LogToConsole':0, 'OutputFlag':0},
               pyomo=False,
               extra_functionality=aux.extra_functionalities,
               store_basis = True
               )
    else:
        n.lopf(solver_name=solver,
               solver_options = {'LogToConsole':0, 'OutputFlag':0},
               pyomo=False,
               extra_functionality=aux.extra_functionalities,
               warmstart = basis_fn,
               store_basis = True
               )
    # ==================================================================================================================
    # Output results
    # ==================================================================================================================
    # output = {
    #     'LCOH (USD/kg)': n.objective/(n.loads.p_set.values[0]/39.4*8760*1000),
    #     'Wind capacity (MW)': n.generators.p_nom_opt['Wind'],
    #     'Solar capacity (MW)': n.generators.p_nom_opt['Solar'],
    #     'Electrolyzer capacity (MW)': n.links.p_nom_opt['Electrolysis'],
    #     'H2 storage capacity (MWh)': n.stores.e_nom_opt['Compressed H2 Store'],
    # }
    # output_array = np.array([
    #     n.objective/(n.loads.p_set.values[0]/39.4*8760*1000),
    #     n.generators.p_nom_opt['Wind'],
    #     n.generators.p_nom_opt['Solar'],
    #     n.links.p_nom_opt['Electrolysis'],
    #     n.stores.e_nom_opt['Compressed H2 Store'],
    # ])
    #!!! for now just return LCOH
    lcoh = n.objective/(n.loads_t.p_set.sum()[0]*39.4*1000) # convert back to kg H2
    print(lcoh)
    basis_fn = n.basis_fn
    return lcoh, basis_fn # need to change this to something that can be added to hexagons

#%% !!! maybe incorporate creating demand profile for each location in for loop
# if trucks will arrive hourly or more to meet quantity of demand, only calculate pipeline
# because trucking costs will be the same. both algorithms cycle through demand center locations
hexagons = gpd.read_file('Data/hexagons_with_country.geojson')

#!!! currently ignoring land-use restrictions-- add later
# excluder = atlite.gis.ExclusionContainer()

cutout = atlite.Cutout("Cutouts/Kenya-2022.nc")

layout = cutout.uniform_layout()
# can add hydro layout here if desired using hydrogen potential map

pv_profile = cutout.pv(
    panel= 'CSi',
    orientation='latitude_optimal',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
pv_profile = pv_profile.rename(dict(dim_0='hexagon'))

wind_profile = cutout.wind(
    turbine = 'Vestas_V80_2MW_gridstreamer',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
wind_profile = wind_profile.rename(dict(dim_0='hexagon'))

investment_excel_path = 'Parameters/investment_parameters.xlsx'
investment_parameters = pd.read_excel(investment_excel_path,
                                    index_col='Country')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  )
demand_centers = demand_parameters.index

for location in demand_centers:
    lcohs_trucking = np.zeros(len(pv_profile.hexagon))
    bases = None
    demand_path=f'Resources/{location} trucking demand.xlsx'
    hydrogen_demand_trucking = pd.read_excel(demand_path,index_col = 0) # Excel file in kg hydrogen, convert to MWh

    start = time.process_time()
    # !!! should probably read demand profile excel here and just input pandas dataframe into
    # function
    for hexagon in pv_profile.hexagon.data:
        investment_series = investment_parameters.loc[hexagons.country[hexagon]]
        if bases == None:
            lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                    pv_profile.sel(hexagon = hexagon),
                                    wind_profile.time,
                                    hydrogen_demand_trucking,
                                    hexagons.loc[hexagon,'theo_turbines'],
                                    hexagons.loc[hexagon,'theo_pv'],
                                    investment_series, 
                                    )
        else:
            # print('Warmstarting...')
            lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                    pv_profile.sel(hexagon = hexagon),
                                    wind_profile.time,
                                    hydrogen_demand_trucking,
                                    hexagons.loc[hexagon,'theo_turbines'],
                                    hexagons.loc[hexagon,'theo_pv'],
                                    investment_series,
                                    basis_fn = bases
                                    )
        lcohs_trucking[hexagon]=lcoh
        bases = basis_fn
    trucking_time = time.process_time()-start
    print(str(trucking_time) + ' s')
    # %% calculate cost of production for pipeline demand profile
    lcohs_pipeline = np.zeros(len(pv_profile.hexagon))
    bases = None
    
    start = time.process_time()
    pipeline_demand_path=f'Resources/{location} pipeline demand.xlsx'
    hydrogen_demand_pipeline = pd.read_excel(pipeline_demand_path, index_col = 0) # Excel file in kg hydrogen, convert to MWh

    for hexagon in pv_profile.hexagon.data:
        investment_series = investment_parameters.loc[hexagons.country[hexagon]]
        if bases == None:
            lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                    pv_profile.sel(hexagon = hexagon),
                                    wind_profile.time,
                                    hydrogen_demand_pipeline,
                                    hexagons.loc[hexagon,'theo_turbines'],
                                    hexagons.loc[hexagon,'theo_pv'],
                                    investment_series,
                                    )
        else:
            # print('Warmstarting...')
            lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                    pv_profile.sel(hexagon = hexagon),
                                    wind_profile.time,
                                    hydrogen_demand_pipeline,
                                    hexagons.loc[hexagon,'theo_turbines'],
                                    hexagons.loc[hexagon,'theo_pv'],
                                    investment_series,
                                    basis_fn = bases
                                    )
        lcohs_pipeline[hexagon]=lcoh
        bases = basis_fn
    pipeline_time = time.process_time()-start
    print(str(pipeline_time) + ' s')

    # add optimal LCOH for each hexagon to hexagon file
    hexagons[f'{location} trucking production cost'] = lcohs_trucking
    hexagons[f'{location} pipeline production cost'] = lcohs_pipeline
    
hexagons.to_file('Resources/hex_lcoh.geojson', driver='GeoJSON')

#%% plot LCOH for each hexagon
# update central coordinates for area considered
crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)

fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'Nairobi trucking production cost',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Production LCOH [euros/kg]'},
    missing_kwds={
        "color": "lightgrey",
        "label": "Missing values",
    },    
)
ax.set_title('Nairobi trucking production cost')


fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'Nairobi pipeline production cost',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Production LCOH [euros/kg]'},
    missing_kwds={
        "color": "lightgrey",
        "label": "Missing values",
    },    
)
ax.set_title('Nairobi pipeline production cost')


fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'Mombasa trucking production cost',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Production LCOH [euros/kg]'},    
    missing_kwds={
        "color": "lightgrey",
        "label": "Missing values",
    },    
)
ax.set_title('Mombasa trucking production cost')
# %% cost difference about 1-4 euro cents per kg H2 or 4-15% more expensive, with a larger difference in the cheapest areas
# at 100 t of annual demand-- difference smaller for larger demands

# hexagons['trucking extra cost'] = hexagons['trucking production cost']/hexagons['pipeline production cost']-1

# fig = plt.figure(figsize=(10,5))

# ax = plt.axes(projection=crs)
# ax.set_axis_off()

# hexagons.to_crs(crs.proj4_init).plot(
#     ax=ax,
#     column = 'trucking extra cost',
#     legend = True,
#     cmap = 'viridis_r',
#     legend_kwds={'label':'Production LCOH [euros/kg]'},    
# )


