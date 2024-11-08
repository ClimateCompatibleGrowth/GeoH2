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
import shapely.wkt
import geopy.distance
import os

def check_folder_exists(name):
    """
    Creates folders if they do not already exist.
    """
    if not os.path.exists(name):
        os.makedirs(name)

def calculate_distance_to_demand(hex_geometry, demand_center_lat, demand_center_lon):
    """
    Calculates distance to demand.

    ...
    Parameters
    ----------
    hex_geometry : geodataframe
        Parameter1 description
    demand_center_lat : 

    demand_center_lon : 

    
    Returns
    -------
    dist : float
        Distance between demand center coordinates and hexagon coordinates in km
    """
    poly = shapely.wkt.loads(str(hex_geometry))
    center = poly.centroid
    demand_coords = (demand_center_lat, demand_center_lon)
    hexagon_coords = (center.y, center.x)
    dist = geopy.distance.geodesic(demand_coords, hexagon_coords).km

    return dist

def calculate_road_construction_cost(distance_to_road, road_capex, 
                                     infrastructure_interest_rate, 
                                     infrastructure_lifetime_years, road_opex):
    """
    Calculates the cost of constructing a road from x to x.

    ...
    Parameters
    ----------
    distance_to_road : float
        
    road_capex : 

    infrastructure_interest_rate : 
    
    infrastructure_lifetime_years :

    road_opex :

    
    Returns
    -------
    cost : float
        Cost to construct a road
    """
    cost = distance_to_road * road_capex * CRF(
                                            infrastructure_interest_rate, 
                                            infrastructure_lifetime_years
                                            ) + distance_to_road * road_opex

    return cost

def main():
    tech_params_filepath = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    conversion_params_filepath = 'parameters/conversion_parameters.xlsx' # SNAKEMAKE INPUT
    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    pipeline_params_filepath = 'parameters/pipeline_parameters.xlsx' # SNAKEMAKE INPUT
    hexagons_filepath = 'data/hex_final_DJ.geojson'


    infra_data = pd.read_excel(tech_params_filepath,
                           sheet_name='Infra',
                           index_col='Infrastructure')
    
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    sheet_name='Demand centers',
                                    index_col='Demand center',
                                    )
    demand_centers = demand_center_list.index
    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')
    

    needs_pipeline_construction = True # CONFIG YAML
    needs_road_construction = True # CONFIG YAML

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']

    hexagons = gpd.read_file(hexagons_filepath)

    check_folder_exists("results")
    
    #%% calculate cost of hydrogen state conversion and transportation for demand
    # loop through all demand centers-- limit this on continential scale
    for demand_center in demand_centers:
        # Demand location based variables
        demand_center_lat = demand_center_list.loc[demand_center,'Lat [deg]']
        demand_center_lon = demand_center_list.loc[demand_center,'Lon [deg]']
        demand_location = Point(demand_center_lon, demand_center_lat)
        hydrogen_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']
        demand_state = demand_center_list.loc[demand_center,'Demand state']

        # Storage hexagons for costs calculated in the next for loop
        road_construction_costs = np.empty(len(hexagons))
        trucking_states = np.empty(len(hexagons),dtype='<U10')
        trucking_costs = np.empty(len(hexagons))
        pipeline_costs = np.empty(len(hexagons))

        # Prices from the country excel file
        elec_price = country_params['Electricity price (euros/kWh)'].iloc[0]
        heat_price = country_params['Heat price (euros/kWh)'].iloc[0]
        plant_interest_rate = country_params['Plant interest rate'].iloc[0]
        infrastructure_interest_rate = country_params['Infrastructure interest rate'].iloc[0]
        infrastructure_lifetime = country_params['Infrastructure lifetime (years)'].iloc[0]
        
        if demand_state not in ['500 bar','LH2','NH3']:
            raise NotImplementedError(f'{demand_state} demand not supported.')
        
        for i in range(len(hexagons)):
            distance_to_road = hexagons['road_dist'][i]
            hex_geometry = hexagons['geometry'][i]
            dist_to_demand = calculate_distance_to_demand(hex_geometry, demand_center_lat, demand_center_lon)

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration
            if hex_geometry.contains(demand_location) == True:
                # Calculate cost of converting hydrogen to a demand state for local demand (i.e. no transport)
                if demand_state == 'NH3':
                    trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state+'_load',
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
                    trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state,
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
            if needs_road_construction == True:
                if distance_to_road==0:
                    road_construction_costs[i] = 0.
                elif distance_to_road!=0 and distance_to_road<10:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(distance_to_road, 
                                                         short_road_capex, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime, 
                                                         road_opex)
                else:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(distance_to_road, 
                                                         long_road_capex, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime, 
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
            elif distance_to_road==0:
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
            elif distance_to_road>0: 
                trucking_costs[i]=trucking_states[i] = np.nan

            # Calculate costs of constructing a pipeline to each hexagon
            if needs_pipeline_construction== True:
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

        # Hexagon file updated with each demand center's costs and states
        hexagons[f'{demand_center} road construction costs'] = road_construction_costs/hydrogen_quantity
        hexagons[f'{demand_center} trucking transport and conversion costs'] = trucking_costs # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} trucking state'] = trucking_states # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} pipeline transport and conversion costs'] = pipeline_costs # cost of supply conversion, pipeline transport, and demand conversion

    hexagons.to_file('results/hex.geojson', driver='GeoJSON', encoding='utf-8') # SNAKEMAKE OUTPUT

if __name__ == "__main__":
    main()

