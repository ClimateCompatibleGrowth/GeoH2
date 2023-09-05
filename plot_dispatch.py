#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 13:57:02 2023

@author: Claire Halloran, University of Oxford

This script plots the temporal results of hydrogen plant optimization for specified
hexagons

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
from optimize_hydrogen_plant import *

logging.basicConfig(level=logging.ERROR)
# list of hexagons of interest
hexagon_list = [372, 9]
# hexagon_list = [9]

# in the future, may want to make hexagons a class with different features
def optimize_hydrogen_plant(wind_potential, pv_potential, times, demand_profile,
                            wind_max_capacity, pv_max_capacity, 
                            country_series, water_limit = None):
    '''
   Optimizes the size of green hydrogen plant components based on renewable potential, hydrogen demand, and country parameters. 

    Parameters
    ----------
    wind_potential : xarray DataArray
        1D dataarray of per-unit wind potential in hexagon.
    pv_potential : xarray DataArray
        1D dataarray of per-unit solar potential in hexagon.
    times : xarray DataArray
        1D dataarray with timestamps for wind and solar potential.
    demand_profile : pandas DataFrame
        hourly dataframe of hydrogen demand in kg.
    country_series : pandas Series
        interest rate and lifetime information.
    water_limit : float
        annual limit on water available for electrolysis in hexagon, in cubic meters. Default is None.
    
    Returns
    -------
    lcoh : float
        levelized cost per kg hydrogen.
    wind_capacity: float
        optimal wind capacity in MW.
    solar_capacity: float
        optimal solar capacity in MW.
    electrolyzer_capacity: float
        optimal electrolyzer capacity in MW.
    battery_capacity: float
        optimal battery storage capacity in MW/MWh (1 hour batteries).
    h2_storage: float
        optimal hydrogen storage capacity in MWh.

    '''    
    
    # if a water limit is given, check if hydrogen demand can be met
    if water_limit != None:
        # total hydrogen demand in kg
        total_hydrogen_demand = demand_profile['Demand'].sum()
        # check if hydrogen demand can be met based on hexagon water availability
        water_constraint =  total_hydrogen_demand <= water_limit * 111.57 # kg H2 per cubic meter of water
        if water_constraint == False:
            print('Not enough water to meet hydrogen demand!')
            # return null values
            lcoh = np.nan
            wind_capacity = np.nan
            solar_capacity = np.nan
            electrolyzer_capacity = np.nan
            battery_capacity = np.nan
            h2_storage = np.nan
            return lcoh, wind_capacity, solar_capacity, electrolyzer_capacity, battery_capacity, h2_storage
    
    # Set up network
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
          p_set = demand_profile['Demand']/1000*39.4,
          )
    
    # Send the weather data to the model
    n.generators_t.p_max_pu['Wind'] = wind_potential
    n.generators_t.p_max_pu['Solar'] = pv_potential
    
    # specify maximum capacity based on land use
    n.generators.loc['Wind','p_nom_max'] = wind_max_capacity
    n.generators.loc['Solar','p_nom_max'] = pv_max_capacity

    # specify technology-specific and country-specific WACC and lifetime here
    n.generators.loc['Wind','capital_cost'] = n.generators.loc['Wind','capital_cost']\
        * CRF(country_series['Wind interest rate'], country_series['Wind lifetime (years)'])
    n.generators.loc['Solar','capital_cost'] = n.generators.loc['Solar','capital_cost']\
        * CRF(country_series['Solar interest rate'], country_series['Solar lifetime (years)'])
    for item in [n.links, n.stores,n.storage_units]:
        item.capital_cost = item.capital_cost * CRF(country_series['Plant interest rate'],country_series['Plant lifetime (years)'])
    
    # Solve the model
    solver = 'gurobi'
    n.lopf(solver_name=solver,
            solver_options = {'LogToConsole':0, 'OutputFlag':0},
            pyomo=False,
            extra_functionality=aux.extra_functionalities,
            )
    # Output results

    lcoh = n.objective/(n.loads_t.p_set.sum()[0]/39.4*1000) # convert back to kg H2
    # wind_capacity = n.generators.p_nom_opt['Wind']
    # solar_capacity = n.generators.p_nom_opt['Solar']
    # electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    # battery_capacity = n.storage_units.p_nom_opt['Battery']
    # h2_storage = n.stores.e_nom_opt['Compressed H2 Store']
    print(lcoh)
    
    
    # instead of saving capacities, want to plot temporal profile of input and output
    # for testing, use time window of July
    return n
    

def plot_dispatch(n, time, demand, hex):
    fig, ax = plt.subplots(figsize=(6, 3))
    
    p_by_carrier = n.generators_t.p.loc[time]
    
    if not n.storage_units.empty:
    
        sto = (
            n.storage_units_t.p.loc[time]
        )
        
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)
    if not n.stores.empty:
        sto = (
            n.stores_t.p.loc[time]
            )
        p_by_carrier = pd.concat([p_by_carrier, sto], axis=1)
    
    p_by_carrier.where(p_by_carrier > 0).loc[time].plot.area(
        ax=ax,
        linewidth=0,
    )
    
    charge = p_by_carrier.where(p_by_carrier < 0).dropna(how="all", axis=1).loc[time]
    
    if not charge.empty:
        charge.plot.area(
            ax=ax,
            linewidth=0,
        )
    
    n.loads_t.p_set.sum(axis=1).loc[time].plot(ax=ax, c="k")
    
    plt.legend(loc=(1.05, 0.2))
    ax.set_ylabel("MW")
    ax.set_ylim(-1500, 2000)
    fig.savefig(f'Resources\\{demand} hexagon {hex} temporal - October.png', bbox_inches='tight')
    plt.close()
    
# %%
time = slice('2022-10-01','2022-10-15')


transport_excel_path = "Parameters/transport_parameters.xlsx"
weather_excel_path = "Parameters/weather_parameters.xlsx"
country_excel_path = 'Parameters/country_parameters.xlsx'
country_parameters = pd.read_excel(country_excel_path,
                                    index_col='Country')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  ).squeeze("columns")
demand_centers = demand_parameters.index
weather_parameters = pd.read_excel(weather_excel_path,
                                   index_col = 'Parameters'
                                   ).squeeze('columns')
weather_filename = weather_parameters['Filename']

hexagons = gpd.read_file('Resources/hex_transport.geojson')
# !!! change to name of cutout in weather
cutout = atlite.Cutout('Cutouts/' + weather_filename +'.nc')
layout = cutout.uniform_layout()

pv_profile = cutout.pv(
    panel= 'CSi',
    orientation='latitude_optimal',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
pv_profile = pv_profile.rename(dict(dim_0='hexagon'))

wind_profile = cutout.wind(
    # Changed turbine type - was Vestas_V80_2MW_gridstreamer in first run
    # Other option being explored: NREL_ReferenceTurbine_2020ATB_4MW, Enercon_E126_7500kW
    turbine = 'NREL_ReferenceTurbine_2020ATB_4MW',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
wind_profile = wind_profile.rename(dict(dim_0='hexagon'))
# %%
for location in demand_centers:
    for hexagon in hexagon_list:
        hydrogen_demand_trucking, hydrogen_demand_pipeline = demand_schedule(
            demand_parameters.loc[location,'Annual demand [kg/a]'],
            hexagons.loc[hexagon,f'{location} trucking state'],
            transport_excel_path,
            weather_excel_path)
        country_series = country_parameters.loc[hexagons.country[hexagon]]
        # trucking demand
        n = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                pv_profile.sel(hexagon = hexagon),
                                wind_profile.time,
                                hydrogen_demand_trucking,
                                hexagons.loc[hexagon,'theo_turbines'],
                                hexagons.loc[hexagon,'theo_pv'],
                                country_series, 
                                # water_limit = hexagons.loc[hexagon,'delta_water_m3']
                                )
        plot_dispatch(n, time, location, hexagon)
        # pipeline demand
        n = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                pv_profile.sel(hexagon = hexagon),
                                wind_profile.time,
                                hydrogen_demand_pipeline,
                                hexagons.loc[hexagon,'theo_turbines'],
                                hexagons.loc[hexagon,'theo_pv'],
                                country_series,
                                # water_limit = hexagons.loc[hexagon,'delta_water_m3'],
                                )
        plot_dispatch(n, time, location, hexagon)

