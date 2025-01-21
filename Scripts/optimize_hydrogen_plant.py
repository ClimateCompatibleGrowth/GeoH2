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
import pandas as pd
import p_H2_aux as aux
from functions import CRF
import numpy as np
import logging
import time

logging.basicConfig(level=logging.ERROR)

def demand_schedule(quantity, start_date, end_date, transport_state, transport_excel_path):
    '''
    calculates hourly hydrogen demand for truck shipment and pipeline transport.

    Parameters
    ----------
    quantity : float
        annual amount of hydrogen to transport in kilograms.
    start_date: string
        start date for demand schedule in the format YYYY-MM-DD.
    end_date: string
        end date for demand schedule in the format YYYY-MM-DD.
    transport_state : string
        state hydrogen is transported in, one of '500 bar', 'LH2', 'LOHC', or 'NH3'.
    transport_excel_path : string
        path to transport_parameters.xlsx file

    Returns
    -------
    trucking_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for hydrogen trucking.
    pipeline_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for pipeline transport.
    '''
    # schedule for pipeline
    index = pd.date_range(start_date, end_date, freq = 'H')
    pipeline_hourly_quantity = quantity/index.size
    pipeline_hourly_demand_schedule = pd.DataFrame(pipeline_hourly_quantity, index=index,  columns = ['Demand'])

    # if demand center is in hexagon
    if transport_state=="None":
        # schedule for trucking
        annual_deliveries = 365*24
        trucking_hourly_demand = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(trucking_hourly_demand, index=index, columns = ['Demand'])
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)

        return trucking_hourly_demand_schedule, pipeline_hourly_demand_schedule
    else:
        transport_parameters = pd.read_excel(transport_excel_path,
                                            sheet_name = transport_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

        truck_capacity = transport_parameters['Net capacity (kg H2)']

        # schedule for trucking
        annual_deliveries = quantity/truck_capacity
        quantity_per_delivery = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(quantity_per_delivery, index=index, columns=['Demand'])
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)

    return trucking_hourly_demand_schedule, pipeline_hourly_demand_schedule

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
    n.add('Load',
          'Hydrogen demand',
          bus = 'Hydrogen',
          p_set = demand_profile['Demand']/1000*39.4,
          )

    # Send the weather data to the model
    n.generators_t.p_max_pu['Wind'] = wind_potential
    n.generators_t.p_max_pu['Solar'] = pv_potential

    # specify maximum capacity based on land use
    n.generators.loc['Wind','p_nom_max'] = wind_max_capacity*4
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
    wind_capacity = n.generators.p_nom_opt['Wind']
    solar_capacity = n.generators.p_nom_opt['Solar']
    electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    battery_capacity = n.storage_units.p_nom_opt['Battery']
    h2_storage = n.stores.e_nom_opt['Compressed H2 Store']
    print(lcoh)
    return lcoh, wind_capacity, solar_capacity, electrolyzer_capacity, battery_capacity, h2_storage


if __name__ == "__main__":
    transport_excel_path = str(snakemake.input.transport_parameters)
    country_excel_path = str(snakemake.input.country_parameters)
    demand_excel_path = str(snakemake.input.demand_parameters)
    country_parameters = pd.read_excel(country_excel_path,
                                        index_col='Country')
    demand_parameters = pd.read_excel(demand_excel_path,
                                      index_col='Demand center',
                                      ).squeeze("columns")
    demand_centers = demand_parameters.index

    weather_year = snakemake.wildcards.weather_year
    end_weather_year = int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    cutout = atlite.Cutout(f'Cutouts/{snakemake.wildcards.country}_{snakemake.wildcards.weather_year}.nc')
    layout = cutout.uniform_layout()

    # can add hydro and other generators here
    pv_profile = cutout.pv(
        panel= 'CSi',
        orientation='latitude_optimal',
        layout = layout,
        shapes = hexagons,
        per_unit = True
        )
    pv_profile = pv_profile.rename(dict(dim_0='hexagon'))
    
    wind_profile = cutout.wind(
        turbine = 'NREL_ReferenceTurbine_2020ATB_4MW',
        layout = layout,
        shapes = hexagons,
        per_unit = True
        )
    wind_profile = wind_profile.rename(dict(dim_0='hexagon'))

    for location in demand_centers:
        # trucking variables
        lcohs_trucking = np.zeros(len(pv_profile.hexagon))
        t_solar_capacities= np.zeros(len(pv_profile.hexagon))
        t_wind_capacities= np.zeros(len(pv_profile.hexagon))
        t_electrolyzer_capacities= np.zeros(len(pv_profile.hexagon))
        t_battery_capacities = np.zeros(len(pv_profile.hexagon))
        t_h2_storages= np.zeros(len(pv_profile.hexagon))
        
        # pipeline variables
        lcohs_pipeline = np.zeros(len(pv_profile.hexagon))
        p_solar_capacities= np.zeros(len(pv_profile.hexagon))
        p_wind_capacities= np.zeros(len(pv_profile.hexagon))
        p_electrolyzer_capacities= np.zeros(len(pv_profile.hexagon))
        p_battery_capacities = np.zeros(len(pv_profile.hexagon))
        p_h2_storages= np.zeros(len(pv_profile.hexagon))

        # function
        for i in pv_profile.hexagon.data:
            hydrogen_demand_trucking, hydrogen_demand_pipeline =\
                demand_schedule(demand_parameters.loc[location,'Annual demand [kg/a]'],
                                start_date,
                                end_date,
                                hexagons.loc[i,f'{location} trucking state'],
                                transport_excel_path)
            
            country_series = country_parameters.loc[hexagons.country[i]]
            
            transport_types = ["trucking", "pipeline"]
            for j in transport_types:
                if j == "trucking":
                    hydrogen_demand = hydrogen_demand_trucking
                else:
                    hydrogen_demand = hydrogen_demand_pipeline

                lcoh, wind_capacity, solar_capacity, electrolyzer_capacity, battery_capacity, h2_storage =\
                    optimize_hydrogen_plant(wind_profile.sel(hexagon = i),
                                            pv_profile.sel(hexagon = i),
                                            wind_profile.time,
                                            hydrogen_demand,
                                            hexagons.loc[i,'theo_turbines'],
                                            hexagons.loc[i,'theo_pv'],
                                            country_series,
                                            # water_limit = hexagons.loc[hexagon,'delta_water_m3']
                                            )
            
                if j == "trucking":
                    lcohs_trucking[i] = lcoh
                    t_solar_capacities[i] = solar_capacity
                    t_wind_capacities[i] = wind_capacity
                    t_electrolyzer_capacities[i] = electrolyzer_capacity
                    t_battery_capacities[i] = battery_capacity
                    t_h2_storages[i] = h2_storage
                else:
                    lcohs_pipeline[i]=lcoh
                    p_solar_capacities[i] = solar_capacity
                    p_wind_capacities[i] = wind_capacity
                    p_electrolyzer_capacities[i] = electrolyzer_capacity
                    p_battery_capacities[i] = battery_capacity
                    p_h2_storages[i] = h2_storage

        # updating trucking hexagons
        hexagons[f'{location} trucking solar capacity'] = t_solar_capacities
        hexagons[f'{location} trucking wind capacity'] = t_wind_capacities
        hexagons[f'{location} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{location} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{location} trucking H2 storage capacity'] = t_h2_storages
        # save trucking LCOH
        hexagons[f'{location} trucking production cost'] = lcohs_trucking            

        # updating pipeline hexagons
        hexagons[f'{location} pipeline solar capacity'] = p_solar_capacities
        hexagons[f'{location} pipeline wind capacity'] = p_wind_capacities
        hexagons[f'{location} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{location} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{location} pipeline H2 storage capacity'] = p_h2_storages

        # add optimal LCOH for each hexagon to hexagon file
        hexagons[f'{location} pipeline production cost'] = lcohs_pipeline

    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')
