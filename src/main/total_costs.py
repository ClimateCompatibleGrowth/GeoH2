#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:44:32 2023

@author: Claire Halloran, University of Oxford

Total hydrogen cost

Bring together all previous data to calculate lowest-cost of commodity
"""

#%% identify lowest-cost strategy: trucking vs. pipeline

import geopandas as gpd
import numpy as np
import pandas as pd
from utils import check_folder_exists

def main():
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    demand_params_filepath = str(snakemake.input.demand_parameters)
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center',
                                    )
    demand_centers = demand_center_list.index
    plant_type = str(snakemake.wildcards.plant_type)

    check_folder_exists("results")

    # Get lowest cost for each transport type
    for demand_center in demand_centers:
        print(f"Calculating total costs for {demand_center} begins...\n")
        if plant_type == "hydrogen":
            trucking_tranport_costs = hexagons[f'{demand_center} trucking transport and conversion costs']
            pipeline_transport_costs = hexagons[f'{demand_center} pipeline transport and conversion costs']
        elif plant_type == "ammonia":
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
            print(f"Calculating total costs for {i+1} of {len(hexagons)}...")
            hexagons.loc[i, f'{demand_center} lowest cost'] = np.nanmin(
                                    [hexagons.loc[i, f'{demand_center} trucking total cost'],
                                    hexagons.loc[i, f'{demand_center} pipeline total cost']
                                    ])
            
    print("\nCalculations complete.\n")
    hexagons.to_file(str(snakemake.output), driver='GeoJSON', encoding='utf-8')

if __name__ == "__main__":
    main()