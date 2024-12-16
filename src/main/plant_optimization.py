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
from network import Network, nh3_pyomo_constraints

def get_demand_schedule(quantity, start_date, end_date, transport_state, transport_params_filepath, freq):
    '''
    Calculates hourly product demand for truck shipment and pipeline transport.

    Parameters
    ----------
    quantity : float
        annual amount of product to transport in kilograms.
    start_date: string
        start date for demand schedule in the format YYYY-MM-DD.
    end_date: string
        end date for demand schedule in the format YYYY-MM-DD.
    transport_state : string
        state product is transported in, one of '500 bar', 'LH2', 'LOHC', or 'NH3'.
    transport_params_filepath : string
        path to transport_parameters.xlsx file

    Returns
    -------
    trucking_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for trucking transport.
    pipeline_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for pipeline transport.
    '''
    # Schedule for pipeline
    # -- Should "freq" be a snakemake input somehow? So the user can define the resolution of the whole run?
    index = pd.date_range(start_date, end_date, freq = freq)
    pipeline_hourly_quantity = quantity/index.size
    pipeline_hourly_demand_schedule = pd.DataFrame(pipeline_hourly_quantity, index=index,  columns = ['Demand'])
    # Resample pipeline schedule
    pipeline_demand_resampled_schedule = pipeline_hourly_demand_schedule.resample(freq).mean()

    # NaN is where there is no road access and no construction so hexagon is infeasible for trucking.
    if pd.isnull(transport_state):
        trucking_hourly_demand_schedule = np.nan
    # If demand center is in hexagon
    elif transport_state=="None":
        # Schedule for trucking
        annual_deliveries = 365*24 # -- This is just the number of hours in a year? Why is this "annual deliveries"?
        # -- Is it just because this is where there is no transport so it creates an even profile?
        # -- I think we need to add a comment to explain anyway.
        trucking_hourly_demand = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(trucking_hourly_demand, index=index, columns = ['Demand'])
        # -- Could we add some comment to explain why this is still called "trucking" schedule despite no trucking happening?
        # First create hourly schedule
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
        # Then resample to desired frequency using mean
        trucking_demand_resampled_schedule = trucking_hourly_demand_schedule.resample(freq).mean()
    else:
        transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = transport_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

        truck_capacity = transport_params['Net capacity (kg of product)']

        # schedule for trucking
        annual_deliveries = quantity/truck_capacity
        quantity_per_delivery = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(quantity_per_delivery, index=index, columns=['Demand'])
        # First create hourly schedule
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
        # Then resample to desired frequency using mean
        trucking_demand_resampled_schedule = trucking_hourly_demand_schedule.resample(freq).mean()

    return trucking_demand_resampled_schedule, pipeline_demand_resampled_schedule

def get_water_constraint(n, demand_profile, water_limit): 
    '''
    Calculates the water constraint.

    Parameters
    ----------
    n : 
        network
    demand_profile : pandas DataFrame
        hourly dataframe of hydrogen demand in kg.
    water_limit : float
        annual limit on water available for electrolysis in hexagon, in cubic meters.

    Returns
    -------
    water_constraint : boolean
        whether there is a water constraint or not.
    '''
    if n.type == "Hydrogen":
        # total hydrogen demand in kg
        total_hydrogen_demand = demand_profile['Demand'].sum()
        # check if hydrogen demand can be met based on hexagon water availability
        water_constraint =  total_hydrogen_demand <= water_limit * 111.57 # kg H2 per cubic meter of water
    elif n.type == "Ammonia":
        # total ammonia demand in kg
        total_ammonia_demand = (
                    (n.loads_t.p_set['Ammonia demand'] * n.snapshot_weightings['objective']).sum() / 6.25 * 1000)
        # total hydrogen demand in kg
        total_hydrogen_demand = total_ammonia_demand * 17 / 3  # convert kg ammonia to kg H2
        # check if hydrogen demand can be met based on hexagon water availability
        water_constraint = total_hydrogen_demand <= water_limit * 111.57  # kg H2 per cubic meter of water
        # note that this constraint is purely stoichiometric-- more water may be needed for cooling or other processes
    return water_constraint

def get_generator_profile(generator, cutout, layout, hexagons, freq):
    '''
    Determines the generation profile of the specified generator in the cutout based on weather data.

    Parameters
    ----------
    generator : string
        The name of the generator type to be used (i.e., "Solar" for pv, "Wind" for wind turbines)
    cutout : 
        A spatial and temporal subset of ERA-5 weather data consisting of grid cells
        https://atlite.readthedocs.io/en/latest/introduction.html
    layout : 
        The capacity to be built in each of the grid_cells.
        https://atlite.readthedocs.io/en/master/ref_api.html
    hexagons :
        Hexagon GeoJSON file.
    
    Returns
    -------
    profile : 
        A profile where weather data has been converted into a generation time-series.
        https://atlite.readthedocs.io/en/master/ref_api.html
    '''
    if generator == "Solar":
        # The panel string should be in the config file as well in case people want to change that in the prep and main.
        # Alycia to double-check that the CSi is 1MW
        profile = cutout.pv(panel= 'CSi',
                            orientation='latitude_optimal',
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True).resample(time=freq).mean()
        profile = profile.rename(dict(dim_0='hexagon'))
    elif generator == "Wind":
        # This string is what we should need to put in the config file (turbine) for both data prep, replacing the constant 4, replacing here.
        profile = cutout.wind(turbine = 'NREL_ReferenceTurbine_2020ATB_4MW',
                            layout = layout,
                            shapes = hexagons,
                            per_unit = True).resample(time=freq).mean()
        profile = profile.rename(dict(dim_0='hexagon'))
    
    return profile

def solve_model(network_class, solver):
    '''
    Solves model using the provided solver.

    Parameters
    ----------
    n : 
        network.
    solver : string
        name of solver to be used.
    '''
    if network_class.type == "Hydrogen":
        network_class.n.lopf(solver_name=solver,
            solver_options = {'LogToConsole':0, 'OutputFlag':0},
            pyomo=False,
            )
    elif network_class.type == "Ammonia":
        network_class.n.lopf(solver_name=solver,
            solver_options={'LogToConsole': 0, 'OutputFlag': 0},
            pyomo=True,
            extra_functionality=nh3_pyomo_constraints,
            )

def get_h2_results(n, generators):
    '''
    Get final results from network optimisation

    Parameters
    ----------
    n : 
        network
    generators : list
        contains types of generators that this plant uses.

    Returns
    -------
    lc : float
        levelized cost per kg product.
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
    lc = n.objective/(n.loads_t.p_set.sum()[0]/39.4*1000) # convert back to kg H2
    for generator in generators:
            generator_capacities[generator] = n.generators.p_nom_opt[f"{generator}"]
    electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    battery_capacity = n.storage_units.p_nom_opt['Battery']
    h2_storage = n.stores.e_nom_opt['Compressed H2 Store']
    return lc, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage

def get_nh3_results(n, generators):
    '''
    Get final results from network optimisation

    Parameters
    ----------
    n : 
        network
    generators : list
        contains types of generators that this plant uses.
    water_limit : float
        annual limit on water available for electrolysis in hexagon, in cubic meters. Default is None.

    Returns
    -------
    lc : float
        levelized cost per kg product.
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
    lc = n.objective / ((n.loads_t.p_set['Ammonia demand'] * n.snapshot_weightings[
        'objective']).sum() / 6.25 * 1000)  # convert back to kg NH3
    for generator in generators:
            generator_capacities[generator] = n.generators.p_nom_opt[generator]
    electrolyzer_capacity = n.links.p_nom_opt['Electrolysis']
    battery_capacity = n.stores.e_nom_opt['Battery']
    h2_storage = n.stores.e_nom_opt['CompressedH2Store']
    # !!! need to save ammonia storage capacity as well
    nh3_storage = n.stores.e_nom_opt['Ammonia']
    return lc, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage, nh3_storage

if __name__ == "__main__":
    # -- Next two lines to be deleted
    # warnings.filterwarnings("ignore")
    logging.basicConfig(level=logging.ERROR)

    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params = pd.read_excel(country_params_filepath, index_col='Country')
    country_series = country_params.iloc[0]
    demand_params = pd.read_excel(demand_params_filepath, index_col='Demand center')
    demand_centers = demand_params.index

    weather_year = 2022             # snakemake.wildcards.weather_year
    end_weather_year = 2023         # int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = "gurobi" # maybe make this into a snakemake wildcard?
    generators = { "Wind" : [], "Solar" : []} # -- this can be gotten rid of whenever config call is added to line 346
    hexagons = gpd.read_file('resources/hex_transport_DJ.geojson') # SNAKEMAKE INPUT
    pipeline_construction = True # snakemake config

    # Get a uniform capacity layout for all grid cells. https://atlite.readthedocs.io/en/master/ref_api.html
    # Alycia to double-check we are using the right layout
    cutout = atlite.Cutout('cutouts/DJ_2022.nc') # SNAKEMAKE INPUT
    layout = cutout.uniform_layout()
    profiles = []
    len_hexagons = len(hexagons)
    water_limit = None # config call
    freq = "3H" # config call
    # freq = "H" # config call
    
    for gen in generators.keys(): # -- use a config call here instead of declaring above
         profiles.append(get_generator_profile(gen, cutout, layout, hexagons, freq))
    
    times = profiles[0].time
    plant_type = "Ammonia" # -- config call
    # plant_type = "Hydrogen" # -- config call

    # Loop through all demand centers -- limit this on continental scale
    for demand_center in demand_centers:
        # Store trucking results
        trucking_lcs = np.zeros(len_hexagons)
        t_generators_capacities = { "Wind" : [], "Solar" : []} # snakemake config
        t_electrolyzer_capacities= np.zeros(len_hexagons)
        t_battery_capacities = np.zeros(len_hexagons)
        t_h2_storages= np.zeros(len_hexagons)
        
        # Store pipeline variables
        pipeline_lcs = np.zeros(len_hexagons)
        p_generators_capacities = { "Wind" : [], "Solar" : []} # snakemake config
        p_electrolyzer_capacities= np.zeros(len_hexagons)
        p_battery_capacities = np.zeros(len_hexagons)
        p_h2_storages= np.zeros(len_hexagons)

        # Ammonia extra storage variables
        if plant_type == "Ammonia":
            t_nh3_storages = np.zeros(len_hexagons)
            p_nh3_storages = np.zeros(len_hexagons)

        annual_demand_quantity = demand_params.loc[demand_center,'Annual demand [kg/a]']

        # Loop through all hexagons
        for i in range(len_hexagons):
            trucking_state = hexagons.loc[i, f'{demand_center} trucking state']
            gen_index = 0
            generators = { "Wind" : [], "Solar" : []} # snakemake config,  already in the config as list, used in map_costs.py line 268
            
            # Get the demand schedule for both pipeline and trucking transport
            trucking_demand_schedule, pipeline_demand_schedule =\
                get_demand_schedule(annual_demand_quantity,
                                start_date,
                                end_date,
                                trucking_state,
                                transport_params_filepath,
                                freq)
            
            # Get the max capacity for each generation type
            for gen in generators.keys():
                if plant_type == "Hydrogen":
                    potential = profiles[gen_index].sel(hexagon = i)
                elif plant_type == "Ammonia":
                    potential = profiles[gen_index].sel(hexagon=i, time=trucking_demand_schedule.index)
                # -- Eventually make a for loop - we can change the theo_turbines name to be Wind
                if gen == "Wind":
                    # -- We'll need to remove this hard-coded 4 eventually CONFIG FILE - 4 MW turbine in spatial data prep
                    max_capacity = hexagons.loc[i,'theo_turbines']*4
                elif gen == "Solar":
                    max_capacity = hexagons.loc[i,'theo_pv']
                # -- Eventually move loops to something like this so we don't have ifs - max_capacity = hexagons.loc[i, gen] * SNAKEMAKE_CONFIG_GEN_SIZE
                
                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                gen_index += 1

            # For each transport type, set up the network and solve
            transport_types = ["trucking", "pipeline"]
            for transport in transport_types:
                network = Network(plant_type, generators) # config call

                # If trucking is viable
                if transport == "trucking" and pd.isnull(trucking_state) == False:
                    network.set_network(trucking_demand_schedule, times, country_series)

                    # Check for water constraint before any solving occurs
                    if water_limit != None:
                        water_constraint = get_water_constraint(network.n, trucking_demand_schedule, water_limit)
                        if water_constraint == False:
                            print('Not enough water to meet demand!')
                            trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                            t_h2_storages[i] = np.nan
                            for gen, capacity in generators_capacities.items():
                                t_generators_capacities[gen].append(np.nan)
                            continue
                    
                    network.set_generators_in_network(country_series)
                    solve_model(network, solver)

                    if plant_type == "Hydrogen": # config call
                        trucking_lcs[i], generators_capacities, \
                        t_electrolyzer_capacities[i], t_battery_capacities[i], \
                        t_h2_storages[i] = get_h2_results(network.n, generators)
                    elif plant_type == "Ammonia": # config call
                        trucking_lcs[i], generators_capacities, \
                        t_electrolyzer_capacities[i], t_battery_capacities[i], \
                        t_h2_storages[i], t_nh3_storages[i]  = get_nh3_results(network.n, generators)
                    
                    for gen, capacity in generators_capacities.items():
                        t_generators_capacities[gen].append(capacity)

                # If the hexagon has no viable trucking state (i.e., no roads reach it), set everything to nan.
                elif transport == "trucking" and pd.isnull(trucking_state) == True:
                    trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                    t_h2_storages[i] = np.nan
                    for gen, capacity in generators_capacities.items():
                        t_generators_capacities[gen].append(np.nan)
                    if plant_type == "Ammonia": # config call 
                        t_nh3_storages[1] = np.nan

                # For pipeline, set it up with pipeline demand schedule if construction is true
                else:
                    if pipeline_construction == True:
                        network.set_network(pipeline_demand_schedule, times, country_series)

                        # Check for water constraint before any solving occurs
                        if water_limit != None:
                            water_constraint = get_water_constraint(network.n, pipeline_demand_schedule, water_limit)
                            if water_constraint == False:
                                print('Not enough water to meet demand!')
                                trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                                t_h2_storages[i] = np.nan
                                for gen, capacity in generators_capacities.items():
                                    t_generators_capacities[gen].append(np.nan)
                                continue

                        network.set_generators_in_network(country_series)
                        solve_model(network, solver)

                        if plant_type == "Hydrogen": # config call
                            pipeline_lcs[i], generators_capacities, \
                            p_electrolyzer_capacities[i], p_battery_capacities[i], \
                            p_h2_storages[i] = get_h2_results(network.n, generators)
                        elif plant_type == "Ammonia": # config call
                            pipeline_lcs[i], generators_capacities, \
                            p_electrolyzer_capacities[i], p_battery_capacities[i], \
                            p_h2_storages[i], p_nh3_storages[i] = get_nh3_results(network.n, generators)

                        for gen, capacity in generators_capacities.items():
                            p_generators_capacities[gen].append(capacity)

                    # If construction is false, you can't transport it - everything gets nan UNLESS you're in the demand centre hexagon
                    else:
                        # -- Demand location has trucking state as None, this is an easy check
                        if trucking_state == "None":
                            network.set_network(pipeline_demand_schedule, times, country_series)

                            # Check for water constraint before any solving occurs
                            if water_limit != None:
                                water_constraint = get_water_constraint(network.n, pipeline_demand_schedule, water_limit)
                                if water_constraint == False:
                                    print('Not enough water to meet demand!')
                                    trucking_lcs[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                                    t_h2_storages[i] = np.nan
                                    for gen, capacity in generators_capacities.items():
                                        t_generators_capacities[gen].append(np.nan)
                                    continue

                            network.set_generators_in_network(country_series)
                            solve_model(network, solver)

                            if plant_type == "Hydrogen": # config call
                                pipeline_lcs[i], generators_capacities, \
                                p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                p_h2_storages[i] = get_h2_results(network.n, generators)
                            elif plant_type == "Ammonia": # config call
                                pipeline_lcs[i], generators_capacities, \
                                p_electrolyzer_capacities[i], p_battery_capacities[i], \
                                p_h2_storages[i], p_nh3_storages[i] = get_nh3_results(network.n, generators)

                            for gen, capacity in generators_capacities.items():
                                p_generators_capacities[gen].append(capacity)
                        else:
                            pipeline_lcs[i], p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = np.nan
                            for gen in p_generators_capacities.keys():
                                p_generators_capacities[gen].append(np.nan)
                            if plant_type == "Ammonia": # config call 
                                t_nh3_storages[1] = np.nan
                
        # Updating trucking-based results in hexagon file
        for gen, capacities in t_generators_capacities.items():
            hexagons[f'{demand_center} trucking {gen.lower()} capacity'] = capacities
        hexagons[f'{demand_center} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{demand_center} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{demand_center} trucking H2 storage capacity'] = t_h2_storages
        hexagons[f'{demand_center} trucking production cost'] = trucking_lcs
        if plant_type == "Ammonia":        
            hexagons[f'{demand_center} pipeline NH3 storage capacity'] = t_nh3_storages

        # Updating pipeline-based results in hexagon file
        for gen, capacities in p_generators_capacities.items():
            hexagons[f'{demand_center} pipeline {gen.lower()} capacity'] = capacities
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{demand_center} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{demand_center} pipeline H2 storage capacity'] = p_h2_storages
        hexagons[f'{demand_center} pipeline production cost'] = pipeline_lcs
        if plant_type == "Ammonia":  
            hexagons[f'{demand_center} pipeline NH3 storage capacity'] = p_nh3_storages


    hexagons.to_file('resources/hex_lc_DJ.geojson', driver='GeoJSON', encoding='utf-8')