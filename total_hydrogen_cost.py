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
import pandas as pd
import numpy as np

hexagons = gpd.read_file('Resources/hex_water.geojson')
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path,
                                  index_col='Demand center',
                                  )

demand_centers = demand_parameters.index
for demand_center in demand_centers:
    hexagons[f'{demand_center} trucking total cost'] =\
        hexagons[f'{demand_center} road construction costs']\
            +hexagons[f'{demand_center} trucking transport and conversion costs']\
                +hexagons[f'{demand_center} trucking production cost']\
                    +hexagons['Lowest water cost']
    hexagons[f'{demand_center} pipeline total cost'] =\
            hexagons[f'{demand_center} pipeline transport and conversion costs']\
                +hexagons[f'{demand_center} pipeline production cost']\
                    +hexagons['Lowest water cost']
                    
    for hexagon in hexagons.index:
        hexagons.loc[hexagon,f'{demand_center} lowest cost'] = np.nanmin(
            [hexagons.loc[hexagon,f'{demand_center} trucking total cost'],
             hexagons.loc[hexagon,f'{demand_center} pipeline total cost']
             ])
        
hexagons.to_file('Resources/hex_total_cost.geojson', driver='GeoJSON', encoding='utf-8')
