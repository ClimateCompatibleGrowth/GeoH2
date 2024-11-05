#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 16:52:59 2023

@author: Claire Halloran, University of Oxford
Contains code originally written by Leander MÃ¼ller, RWTH Aachen University

Calculates the cost-optimal hydrogen transportation strategy to the nearest demand center.

Calculate cost of pipeline transport and demand profile based on optimal size



"""

import geopandas as gpd
import numpy as np
import pandas as pd
from functions import CRF, cheapest_trucking_strategy, h2_conversion_stand, cheapest_pipeline_strategy
from shapely.geometry import Point
import shapely.geometry
import shapely.wkt
import geopy.distance
import os
import json

def check_folder_exists(name):
    """
    Create folders if they do not already exist.
    """
    if not os.path.exists(name):
        os.makedirs(name)

def calculate_distance_to_demand(hexagon, i, demand_center_list, d):
    """
    Calculate distance to demand.

    ...
    Parameters
    ----------
    hexagon : geodataframe
        Parameter1 description
    i : int
        Parameter2 description
    demand_center_list : 

    d : 

    
    Returns
    -------
    dist : float
        Distance between demand center coordinates and hexagon coordinates in km
    """
    poly = shapely.wkt.loads(str(hexagon['geometry'][i]))
    center = poly.centroid
    demand_coords = (demand_center_list.loc[d,'Lat [deg]'], demand_center_list.loc[d,'Lon [deg]'])
    hexagon_coords = (center.y, center.x)
    dist = geopy.distance.geodesic(demand_coords, hexagon_coords).km

    return dist

def calculate_road_construction_cost(road_distance, road_capex, 
                                     infrastructure_interest_rate, 
                                     infrastructure_lifetime_years, road_opex):
    """
    Calculate road construction cost.

    ...
    Parameters
    ----------
    road_distance : float
        
    road_capex : 

    infrastructure_interest_rate : 
    
    infrastructure_lifetime_years :

    road_opex :

    
    Returns
    -------
    cost : float
        Cost to construct a road
    """
    cost = road_distance * road_capex * CRF(
                                            infrastructure_interest_rate, 
                                            infrastructure_lifetime_years
                                            ) + road_distance * road_opex

    return cost

def main():
    tech_params_filepath = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    conversion_params_filepath = 'parameters/conversion_parameters.xlsx' # SNAKEMAKE INPUT
    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    pipeline_params_filepath = 'parameters/pipeline_parameters.xlsx' # SNAKEMAKE INPUT
    hexagon_filepath = 'data/hex_final_DJ.geojson'


    infra_data = pd.read_excel(tech_params_filepath,
                           sheet_name='Infra',
                           index_col='Infrastructure')
    
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    sheet_name='Demand centers',
                                    index_col='Demand center',
                                    )
    country_parameters = pd.read_excel(country_params_filepath,
                                        index_col='Country')

    pipeline_construction = True # CONFIG YAML
    road_construction = True # CONFIG YAML

    road_capex_long = infra_data.at['Long road','CAPEX']
    road_capex_short = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']

    hexagon = gpd.read_file(hexagon_filepath)

    check_folder_exists("resources")
    
    #%% calculate cost of hydrogen state conversion and transportation for demand
    # loop through all demand centers-- limit this on continential scale
    for d in demand_center_list.index:
        # Demand location based variables
        demand_location = Point(demand_center_list.loc[d,'Lon [deg]'], demand_center_list.loc[d,'Lat [deg]'])
        hydrogen_quantity = demand_center_list.loc[d,'Annual demand [kg/a]']
        demand_state = demand_center_list.loc[d,'Demand state']

        # Storage hexagons for costs calculated in the next for loop
        road_construction_costs = np.empty(len(hexagon))
        trucking_states = np.empty(len(hexagon),dtype='<U10')
        trucking_costs = np.empty(len(hexagon))
        pipeline_costs = np.empty(len(hexagon))

        # Prices from the country excel file
        elec_price = country_parameters['Electricity price (euros/kWh)'].iloc[0]
        heat_price = country_parameters['Heat price (euros/kWh)'].iloc[0]
        plant_interest_rate = country_parameters['Plant interest rate'].iloc[0]
        infrastructure_interest_rate = country_parameters['Infrastructure interest rate'].iloc[0]
        infrastructure_lifetime_years = country_parameters['Infrastructure lifetime (years)'].iloc[0]
        
        if demand_state not in ['500 bar','LH2','NH3']:
            raise NotImplementedError(f'{demand_state} demand not supported.')
        
        for i in range(len(hexagon)):
            dist_to_demand = calculate_distance_to_demand(hexagon, i, demand_center_list, d)
            road_distance = hexagon['road_dist'][i]

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration
            if hexagon['geometry'][i].contains(demand_location) == True:
                # Calculate cost of converting hydrogen to a demand state for local demand (i.e. no transport)
                if demand_state == 'NH3':
                    trucking_costs[i], pipeline_costs[i] = h2_conversion_stand(demand_state+'_load',
                                             hydrogen_quantity,
                                             elec_price,
                                             heat_price,
                                             plant_interest_rate,
                                             conversion_params_filepath
                                             )[2]/hydrogen_quantity
                    trucking_states[i] = "None"
                    road_construction_costs[i] = 0.
                    continue
                else:
                    trucking_costs[i], pipeline_costs[i] = h2_conversion_stand(demand_state,
                                             hydrogen_quantity,
                                             elec_price,
                                             heat_price,
                                             plant_interest_rate,
                                             conversion_params_filepath
                                             )[2]/hydrogen_quantity
                    trucking_states[i] = "None"
                    road_construction_costs[i] = 0.
                    continue

            # determine elec_cost at demand to determine potential energy costs
            # Calculate cost of constructing a road to each hexagon
            if road_construction == True:
                print(hexagon['road_dist'][i])
                if road_distance==0:
                    road_construction_costs[i] = 0.
                elif road_distance!=0 and road_distance<10:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(road_distance, 
                                                         road_capex_short, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime_years, 
                                                         road_opex)
                else:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(road_distance, 
                                                         road_capex_long, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime_years, 
                                                         road_opex)
                    
                trucking_costs[i], trucking_states[i] =\
                    cheapest_trucking_strategy(demand_state,
                                                hydrogen_quantity,
                                                dist_to_demand,
                                                elec_price,
                                                heat_price,
                                                infrastructure_interest_rate,
                                                conversion_params_filepath,
                                                transport_params_filepath,
                                                )
            elif road_distance==0:
                trucking_costs[i], trucking_states[i] =\
                    cheapest_trucking_strategy(demand_state,
                                                hydrogen_quantity,
                                                dist_to_demand,
                                                elec_price,
                                                heat_price,
                                                infrastructure_interest_rate,
                                                conversion_params_filepath,
                                                transport_params_filepath,
                                                )
            elif road_distance>0: 
                trucking_costs[i] = np.nan
                trucking_states[i] = np.nan

            # Calculate costs of constructing a pipeline to each hexagon
            if pipeline_construction== True:
                pipeline_cost, pipeline_type =\
                    cheapest_pipeline_strategy(demand_state,
                                            hydrogen_quantity,
                                            dist_to_demand,
                                            elec_price,
                                            heat_price,
                                            infrastructure_interest_rate,
                                            conversion_params_filepath,
                                            pipeline_params_filepath,
                                            )
                pipeline_costs[i] = pipeline_cost
            else:
                pipeline_costs[i] = np.nan

        # variables to save for each demand scenario
        hexagon[f'{d} road construction costs'] = road_construction_costs/hydrogen_quantity
        hexagon[f'{d} trucking transport and conversion costs'] = trucking_costs # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagon[f'{d} trucking state'] = trucking_states # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagon[f'{d} pipeline transport and conversion costs'] = pipeline_costs # cost of supply conversion, pipeline transport, and demand conversion

    # Added force to UTF-8 encoding.
    hexagon.to_file('resources/hex_transport_DJ.geojson', driver='GeoJSON', encoding='utf-8') # SNAKEMAKE OUTPUT

if __name__ == "__main__":
    main()

