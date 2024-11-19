# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:47:57 2023

@author: Claire Halloran, University of Oxford

Includes code from Nicholas Salmon, University of Oxford, for optimizing
hydrogen plant capacity.

"""

import atlite
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import warnings
from network import Network

def get_demand_schedule(quantity, start_date, end_date, transport_state, transport_params_filepath):
    '''
    Calculates hourly hydrogen demand for truck shipment and pipeline transport.

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
    transport_params_filepath : string
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
    else:
        transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = transport_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

        truck_capacity = transport_params['Net capacity (kg H2)']

        # schedule for trucking
        annual_deliveries = quantity/truck_capacity
        quantity_per_delivery = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(quantity_per_delivery, index=index, columns=['Demand'])
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)

    return trucking_hourly_demand_schedule, pipeline_hourly_demand_schedule

def get_generator_profile(generator, cutout, layout, hexagons):
    '''
    Sets the profile of the specified generator in the cutout.

    Parameters
    ----------
    generator : string
        ...
    cutout : 
        ...
    layout : 
        ...
    hexagons :
        ...
    Returns
    -------
    profile : 
        ...
    '''
    if generator == "Solar":
        profile = cutout.pv(panel= 'CSi',
                            orientation='latitude_optimal',
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True)
        profile = profile.rename(dict(dim_0='hexagon'))
    elif generator == "Wind":
        profile = cutout.wind(turbine = 'NREL_ReferenceTurbine_2020ATB_4MW',
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True)
        profile = profile.rename(dict(dim_0='hexagon'))
    
    return profile

def solve_model(n, solver):
    '''
    Solves model using the provided solver.

    Parameters
    ----------
    n : 
        network
    solver : string
        name of solver to be used.
    '''
    n.lopf(solver_name=solver,
        solver_options = {'LogToConsole':0, 'OutputFlag':0},
        pyomo=False,
        )

def _get_water_constraint(demand_profile, water_limit): 
    '''
    Calculates the water constraint.

    Parameters
    ----------
    demand_profile : pandas DataFrame
        hourly dataframe of hydrogen demand in kg.
    water_limit : float
        annual limit on water available for electrolysis in hexagon, in cubic meters.

    Returns
    -------
    water_constraint : boolean
        whether there is a water constraint or not.
    '''
    # total hydrogen demand in kg
    total_hydrogen_demand = demand_profile['Demand'].sum()
    # check if hydrogen demand can be met based on hexagon water availability
    water_constraint =  total_hydrogen_demand <= water_limit * 111.57 # kg H2 per cubic meter of water
    
    return water_constraint

def get_results(n, demand_profile, generators, water_limit=None):
    '''
    Calculates the water constraint.

    Parameters
    ----------
    demand_profile : pandas DataFrame
        hourly dataframe of hydrogen demand in kg.
    generators : list
        contains types of generators that this plant uses.

    Returns
    -------
    lcoh : float
        levelized cost per kg hydrogen.
    generator_capactities : dictionary
        contains each generator with their optimal capacity in MW.
    electrolyzer_capacity: float
        optimal electrolyzer capacity in MW.
    battery_capacity: float
        optimal battery storage capacity in MW/MWh (1 hour batteries).
    h2_storage: float
        optimal hydrogen storage capacity in MWh.
    '''
    generator_capacities = {}
    if water_limit != None :
        water_constraint = _get_water_constraint(demand_profile, water_limit)
        if water_constraint == False:
                    print('Not enough water to meet hydrogen demand!')
                    lcoh = np.nan
                    for generator in generators:
                            generator_capacities[generator] = np.nan
                    electrolyzer_capacity = np.nan
                    battery_capacity = np.nan
                    h2_storage = np.nan

    if water_limit == None :
        lcoh = n.objective/(n.loads_t.p_set.sum()[0]/39.4*1000) # convert back to kg H2
        for generator in generators:
                generator_capacities[generator] = n.generators.p_nom_opt[f"{generator}"]
        electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
        battery_capacity = n.storage_units.p_nom_opt['Battery']
        h2_storage = n.stores.e_nom_opt['Compressed H2 Store']

    return lcoh, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    logging.basicConfig(level=logging.ERROR)

    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')
    country_series = country_params.iloc[0]
    demand_params = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center',
                                    )
    demand_centers = demand_params.index

    weather_year = 2022             # snakemake.wildcards.weather_year
    end_weather_year = 2023         # int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = "gurobi" # maybe make this into a snakemake wildcard?
    generators = { "Wind" : [], "Solar" : []} # already in the config as list, used in map_costs.py line 268
    hexagons = gpd.read_file('results/completed_hex_DJ.geojson') # SNAKEMAKE INPUT

    cutout = atlite.Cutout('cutouts/DJ_2022.nc') # SNAKEMAKE INPUT
    layout = cutout.uniform_layout()
    profiles = []
    
    for gen in generators.keys():
         profiles.append(get_generator_profile(gen, cutout, layout, hexagons))
    
    times = profiles[0].time
    # pv_profile = get_generator_profile("Solar", cutout, layout, hexagons)
    # wind_profile = get_generator_profile("Wind", cutout, layout, hexagons)

    for demand_center in demand_centers:
        # trucking variables
        lcohs_trucking = np.zeros(len(hexagons))
        t_solar_capacities= np.zeros(len(hexagons))
        t_wind_capacities= np.zeros(len(hexagons))
        t_electrolyzer_capacities= np.zeros(len(hexagons))
        t_battery_capacities = np.zeros(len(hexagons))
        t_h2_storages= np.zeros(len(hexagons))
        
        # pipeline variables
        lcohs_pipeline = np.zeros(len(hexagons))
        p_solar_capacities= np.zeros(len(hexagons))
        p_wind_capacities= np.zeros(len(hexagons))
        p_electrolyzer_capacities= np.zeros(len(hexagons))
        p_battery_capacities = np.zeros(len(hexagons))
        p_h2_storages= np.zeros(len(hexagons))

        hydrogen_quantity = demand_params.loc[demand_center,'Annual demand [kg/a]']

        for i in range(len(hexagons)):
            trucking_state = hexagons.loc[i, f'{demand_center} trucking state']
            gen_index = 0
            if i > 0:
                     generators = { "Wind" : [], "Solar" : []}
            for gen in generators.keys():
                potential = profiles[gen_index].sel(hexagon = i)
                if gen == "Wind":
                    max_capacity = hexagons.loc[i,'theo_turbines']*4
                elif gen == "Solar":
                    max_capacity = hexagons.loc[i,'theo_pv']

                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                gen_index += 1
            
            # print(generators)
            # wind_potential = wind_profile.sel(hexagon = i)
            # pv_potential = pv_profile.sel(hexagon = i)
            # wind_max_capacity = hexagons.loc[i,'theo_turbines']*4
            # pv_max_capacity = hexagons.loc[i,'theo_pv']
            # generators = {
            #         "Wind" : [wind_potential, wind_max_capacity],
            #         "Solar" : [pv_potential, pv_max_capacity]
            #      }
            
            trucking_demand_schedule, pipeline_demand_schedule =\
                get_demand_schedule(hydrogen_quantity,
                                start_date,
                                end_date,
                                trucking_state,
                                transport_params_filepath)
            
            transport_types = ["trucking", "pipeline"]
            for transport in transport_types:
                network = Network("Hydrogen", generators)
                if transport == "trucking":
                    network.set_network(trucking_demand_schedule, times, country_series)
                    network.set_generators_in_network(country_series)
                    solve_model(network.n, solver)
                    lcohs_trucking[i], generator_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], t_h2_storages[i] = get_results(network.n, trucking_demand_schedule, generators)
                    t_solar_capacities[i] = generator_capacities["Solar"]
                    t_wind_capacities[i] = generator_capacities["Wind"]
                else:
                    network.set_network(pipeline_demand_schedule, times, country_series)
                    network.set_generators_in_network(country_series)
                    solve_model(network.n, solver)
                    lcohs_pipeline[i], generator_capacities, p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = get_results(network.n, pipeline_demand_schedule, generators)
                    p_solar_capacities[i] = generator_capacities["Solar"]
                    p_wind_capacities[i] = generator_capacities["Wind"]
                
        # updating trucking hexagons
        # LOOK HERE
        for gen in generators:
             pass
        hexagons[f'{demand_center} trucking solar capacity'] = t_solar_capacities
        hexagons[f'{demand_center} trucking wind capacity'] = t_wind_capacities
        hexagons[f'{demand_center} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{demand_center} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{demand_center} trucking H2 storage capacity'] = t_h2_storages
        # save trucking LCOH
        hexagons[f'{demand_center} trucking production cost'] = lcohs_trucking            

        # updating pipeline hexagons
        hexagons[f'{demand_center} pipeline solar capacity'] = p_solar_capacities
        hexagons[f'{demand_center} pipeline wind capacity'] = p_wind_capacities
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{demand_center} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{demand_center} pipeline H2 storage capacity'] = p_h2_storages
        # add optimal LCOH for each hexagon to hexagon file
        hexagons[f'{demand_center} pipeline production cost'] = lcohs_pipeline


    hexagons.to_file("results/hex.geojson", driver='GeoJSON', encoding='utf-8')