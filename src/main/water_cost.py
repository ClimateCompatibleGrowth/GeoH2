#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:26:19 2023

@author: Claire Halloran, University of Oxford

Water costs for hydrogen production in each hexagon


"""

import geopandas as gpd
import pandas as pd
import numpy as np

def main():
    hexagons = gpd.read_file('Resources/hex_transport_DJ.geojson') # SNAKEMAKE INPUT
    tech_params_filepath = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT

    water_data = pd.read_excel(tech_params_filepath,
                                sheet_name='Water',
                                index_col='Parameter'
                                ).squeeze("columns")
    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')

    #%% water cost for each hexagon for each kg hydrogen produced

    # Water cost related variables
    h2o_costs_dom_water_bodies = np.empty(len(hexagons))
    h2o_costs_ocean = np.empty(len(hexagons))
    min_h2o_costs = np.empty(len(hexagons))

    electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
    electricity_demand_ocean_h2o_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
    water_transport_costs = water_data['Water transport cost (euros/100 km/m3)']
    water_spec_cost = water_data['Water specific cost (euros/m3)']
    water_demand = water_data['Water demand  (L/kg H2)']
    elec_price = country_params['Electricity price (euros/kWh)'].iloc[0]
    
    # Calculating water costs for each hexagon
    for i in range(len(hexagons)):
        waterbody_dist = hexagons['waterbody_dist'][i]
        waterway_dist = hexagons['waterway_dist'][i]
        ocean_dist = hexagons['ocean_dist'][i]
        
        h2o_costs_dom_water_bodies[i] =(water_spec_cost 
                                            + (water_transport_costs/100)
                                            * min(waterbody_dist, waterway_dist) 
                                            + electricity_demand_h2o_treatment
                                            * elec_price
                                            ) * water_demand/1000
        h2o_costs_ocean[i] =(water_spec_cost 
                                + (water_transport_costs/100)
                                * ocean_dist
                                + electricity_demand_ocean_h2o_treatment
                                * elec_price
                                ) * water_demand/1000
        
        min_h2o_costs[i] = min(h2o_costs_dom_water_bodies[i], h2o_costs_ocean[i])

    hexagons['Ocean water costs'] = h2o_costs_ocean
    hexagons['Freshwater costs'] = h2o_costs_dom_water_bodies
    hexagons['Lowest water cost'] = min_h2o_costs

    hexagons.to_file("Resources/hex_transport_DJ.geojson", driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()