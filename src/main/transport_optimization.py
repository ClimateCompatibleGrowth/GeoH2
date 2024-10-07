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
        Distance between x and x
    """
    poly = shapely.wkt.loads(str(hexagon['geometry'][i]))
    center = poly.centroid
    demand_coords = (demand_center_list.loc[d,'Lat [deg]'], demand_center_list.loc[d,'Lon [deg]'])
    hexagon_coords = (center.y, center.x)
    dist = geopy.distance.geodesic(demand_coords, hexagon_coords).km

    return dist

def main():
    # Paths
    technology_parameters_path = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    demand_parameters_path = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_parameters_path = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    conversion_parameters_path = 'parameters/conversion_parameters.xlsx' # SNAKEMAKE INPUT
    transport_parameters_path = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    pipeline_parameters_path = 'parameters/pipeline_parameters.xlsx' # SNAKEMAKE INPUT
    hexagon_path = 'data/hex_final_DJ.geojson'

    # Input parameters
    infra_data = pd.read_excel(technology_parameters_path,
                           sheet_name='Infra',
                           index_col='Infrastructure')

    demand_center_list = pd.read_excel(demand_parameters_path,
                                    sheet_name='Demand centers',
                                    index_col='Demand center',
                                    )
    country_parameters = pd.read_excel(country_parameters_path,
                                        index_col='Country')

    pipeline_construction = True # CONFIG YAML
    road_construction = True # CONFIG YAML

    road_capex_long = infra_data.at['Long road','CAPEX']
    road_capex_short = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']

    hexagon = gpd.read_file(hexagon_path)

    check_folder_exists("Resources")

    #%% calculate cost of hydrogen state conversion and transportation for demand
    # loop through all demand centers-- limit this on continential scale
    for d in demand_center_list.index:
        demand_location = Point(demand_center_list.loc[d,'Lon [deg]'], demand_center_list.loc[d,'Lat [deg]'])
        # distance_to_demand = np.empty(len(hexagon))
        hydrogen_quantity = demand_center_list.loc[d,'Annual demand [kg/a]']

        # road_construction_costs = np.empty(len(hexagon))
        # trucking_states = np.empty(len(hexagon),dtype='<U10')
        # trucking_costs = np.empty(len(hexagon))
        # pipeline_costs = np.empty(len(hexagon))

        demand_state = demand_center_list.loc[d,'Demand state']
        demand_fid = 0
        if demand_state not in ['500 bar','LH2','NH3']:
            raise NotImplementedError(f'{demand_state} demand not supported.')

        for i in range(len(hexagon)):
            dist = calculate_distance_to_demand(hexagon, i, demand_center_list, d)
            
            # distance_to_demand[i] = dist

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration
            if hexagon['geometry'][i].contains(demand_location) == True:
                ## label demand location under consideration
                #demand_fid = i

                # calculate cost of converting hydrogen to ammonia for local demand (i.e. no transport)
                if demand_state == 'NH3':
                    local_conversion_cost =\
                        h2_conversion_stand(demand_state+'_load',
                                            hydrogen_quantity,
                                            country_parameters.loc[hexagon['country'][i], 'Electricity price (euros/kWh)'],
                                            country_parameters.loc[hexagon['country'][i], 'Heat price (euros/kWh)'],
                                            country_parameters.loc[hexagon['country'][i],'Plant interest rate'],
                                            conversion_parameters_path
                                            )[2]/hydrogen_quantity

                    hexagon[f'{d} trucking transport and conversion costs'] = local_conversion_cost
                    hexagon[f'{d} pipeline transport and conversion costs'] = local_conversion_cost
                    hexagon[f'{d} trucking state'] = "None"
                    hexagon[f'{d} road construction costs'] = 0.
                    continue
                else:
                    local_conversion_cost =\
                        h2_conversion_stand(demand_state,
                                            hydrogen_quantity,
                                            country_parameters.loc[hexagon['country'][i], 'Electricity price (euros/kWh)'],
                                            country_parameters.loc[hexagon['country'][i], 'Heat price (euros/kWh)'],
                                            country_parameters.loc[hexagon['country'][i],'Plant interest rate'],
                                            conversion_parameters_path
                                            )[2]/hydrogen_quantity
                    hexagon[f'{d} trucking transport and conversion costs'] = local_conversion_cost
                    hexagon[f'{d} pipeline transport and conversion costs'] = local_conversion_cost
                    hexagon[f'{d} trucking state'] = "None"
                    hexagon[f'{d} road construction costs'] = 0.
                    continue

            # determine elec_cost at demand to determine potential energy costs
            # calculate cost of constructing a road to each hexagon
            if road_construction == True:
                if hexagon['road_dist'][i]==0:
                    hexagon[f'{d} road construction costs'] = 0.
                elif hexagon['road_dist'][i]!=0 and hexagon['road_dist'][i]<10:
                    hexagon[f'{d} road construction costs'] = hexagon['road_dist'][i]\
                        *road_capex_short*CRF(
                            country_parameters.loc[hexagon['country'][i],'Infrastructure interest rate'],
                            country_parameters.loc[hexagon['country'][i],'Infrastructure lifetime (years)'])\
                        +hexagon['road_dist'][i]*road_opex
                else:
                    hexagon[f'{d} road construction costs'] = hexagon['road_dist'][i]*road_capex_long*CRF(
                        country_parameters.loc[hexagon['country'][i],'Infrastructure interest rate'],
                        country_parameters.loc[hexagon['country'][i],'Infrastructure lifetime (years)'])\
                    +hexagon['road_dist'][i]*road_opex
                    
                trucking_cost, trucking_state =\
                    cheapest_trucking_strategy(demand_state,
                                                hydrogen_quantity,
                                                dist,
                                                country_parameters.loc[hexagon.country[i],'Electricity price (euros/kWh)'],
                                                country_parameters.loc[hexagon.country[i],'Heat price (euros/kWh)'],
                                                country_parameters.loc[hexagon['country'][i],'Infrastructure interest rate'],
                                                conversion_parameters_path,
                                                transport_parameters_path,
                                                country_parameters.loc[hexagon.country[demand_fid],'Electricity price (euros/kWh)'],
                                                )
                hexagon[f'{d} trucking transport and conversion costs'] = trucking_cost
                hexagon[f'{d} trucking state'] = trucking_state

            elif hexagon['road_dist'][i]==0:
                trucking_cost, trucking_state =\
                    cheapest_trucking_strategy(demand_state,
                                                hydrogen_quantity,
                                                dist,
                                                country_parameters.loc[hexagon.country[i],'Electricity price (euros/kWh)'],
                                                country_parameters.loc[hexagon.country[i],'Heat price (euros/kWh)'],
                                                country_parameters.loc[hexagon['country'][i],'Infrastructure interest rate'],
                                                conversion_parameters_path,
                                                transport_parameters_path,
                                                country_parameters.loc[hexagon.country[demand_fid],'Electricity price (euros/kWh)'],
                                                )
                hexagon[f'{d} trucking transport and conversion costs'] = trucking_cost
                hexagon[f'{d} trucking state'] = trucking_state

            elif hexagon['road_dist'][i]>0: 
                hexagon[f'{d} trucking transport and conversion costs'] = np.nan
                hexagon[f'{d} trucking state'] = np.nan
            # pipeline costs
            if pipeline_construction== True:
            
                pipeline_cost, pipeline_type =\
                    cheapest_pipeline_strategy(demand_state,
                                            hydrogen_quantity,
                                            dist,
                                            country_parameters.loc[hexagon.country[i],'Electricity price (euros/kWh)'],
                                            country_parameters.loc[hexagon.country[i],'Heat price (euros/kWh)'],
                                            country_parameters.loc[hexagon['country'][i],'Infrastructure interest rate'],
                                            conversion_parameters_path,
                                            pipeline_parameters_path,
                                            country_parameters.loc[hexagon.country[demand_fid],'Electricity price (euros/kWh)'],
                                            )
                hexagon[f'{d} pipeline transport and conversion costs'] = pipeline_cost
            else:
                hexagon[f'{d} pipeline transport and conversion costs'] = np.nan

        # variables to save for each demand scenario
        # hexagon[f'{d} road construction costs'] = road_construction_costs/hydrogen_quantity
        # hexagon[f'{d} trucking transport and conversion costs'] = trucking_costs # cost of road construction, supply conversion, trucking transport, and demand conversion
        # hexagon[f'{d} trucking state'] = trucking_states # cost of road construction, supply conversion, trucking transport, and demand conversion
        # hexagon[f'{d} pipeline transport and conversion costs'] = pipeline_costs # cost of supply conversion, pipeline transport, and demand conversion

    # Added force to UTF-8 encoding.
    hexagon.to_file('Resources/hex_transport_DJ.geojson', driver='GeoJSON', encoding='utf-8') # SNAKEMAKE OUTPUT

if __name__ == "__main__":
    main()

