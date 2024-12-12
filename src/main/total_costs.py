#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:44:32 2023

@author: Claire Halloran, University of Oxford

Total hydrogen cost

Bring together all previous data to calculate lowest-cost hydrogen
"""

#%% identify lowest-cost strategy: trucking vs. pipeline

import geopandas as gpd
import numpy as np
import pandas as pd
from utils import check_folder_exists

def main():
    hexagons = gpd.read_file('resources/hex_water_DJ.geojson') # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center',
                                    )
    demand_centers = demand_center_list.index
    plant_type = "Ammonia" # config file

    check_folder_exists("results")

    # Get lowest cost for each transport type
    for demand_center in demand_centers:
        if plant_type == "Hydrogen":
            trucking_tranport_costs = hexagons[f'{demand_center} trucking transport and conversion costs']
            pipeline_transport_costs = hexagons[f'{demand_center} pipeline transport and conversion costs']
        elif plant_type == "Ammonia":
            trucking_tranport_costs = hexagons[f'{demand_center} trucking transport costs']
            pipeline_transport_costs = hexagons[f'{demand_center} pipeline transport costs']

        hexagons[f'{demand_center} trucking total cost'] =\
            hexagons[f'{demand_center} road construction costs'] +\
                trucking_tranport_costs +\
                    hexagons[f'{demand_center} trucking production cost'] +\
                        hexagons['Lowest water cost']
        hexagons[f'{demand_center} pipeline total cost'] =\
                pipeline_transport_costs +\
                    hexagons[f'{demand_center} pipeline production cost'] +\
                        hexagons['Lowest water cost']

        # Get the lowest between the trucking and the pipeline options
        for i in range(len(hexagons)):
            hexagons.loc[i, f'{demand_center} lowest cost'] = np.nanmin(
                                    [hexagons.loc[i, f'{demand_center} trucking total cost'],
                                    hexagons.loc[i, f'{demand_center} pipeline total cost']
                                    ])
            
    hexagons.to_file("results/hex_total_cost_DJ_2022.geojson", driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()