"""
Created on Tue Aug 29th 2023

@author: Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

Cost by hydrogen plant component

Add attributes to hex file for cost of each component
"""

# %% identify lowest-cost strategy: trucking vs. pipeline

import geopandas as gpd
import pandas as pd
from geopy.geocoders import Photon
import functions

# Load hexagons
hexagons = gpd.read_file(str(snakemake.input.hexagons))

# Load necessary parameters
demand_excel_path = str(snakemake.input.demand_parameters)
demand_parameters = pd.read_excel(demand_excel_path, index_col='Demand center')
country_excel_path = str(snakemake.input.country_parameters)
country_parameters = pd.read_excel(country_excel_path, index_col='Country')
stores_csv_path = str(snakemake.input.stores_parameters) # H2 storage
stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
storage_csv_path = str(snakemake.input.storage_parameters) # Battery
storage_parameters = pd.read_csv(storage_csv_path, index_col='name')
links_csv_path = str(snakemake.input.links_parameters) # Electrolyzer
links_parameters = pd.read_csv(links_csv_path, index_col='name')
generators_csv_path = str(snakemake.input.generators_parameters) # Solar and generator
generators_parameters = pd.read_csv(generators_csv_path, index_col='name')

# For each demand center, get costs for each component

demand_centers = demand_parameters.index
transport_methods = ['pipeline', 'trucking']
for demand_center in demand_centers:
    # Get location of demand center
    lat = demand_parameters.loc[demand_center, 'Lat [deg]']
    lon = demand_parameters.loc[demand_center, 'Lon [deg]']
    coordinates = str(lat) + ", " + str(lon)
    # Get country where the demand center is
    geolocator = Photon(user_agent="MyApp")
    location = geolocator.reverse(coordinates, language="en")
    country = location.raw['properties']['country']
    
    # Get CRF and then cost for each component using the data for the country you are looking at
    for transport_method in transport_methods:
        # Battery 
        interest_battery = country_parameters.loc[country, 'Plant interest rate']
        lifetime_battery = country_parameters.loc[country, 'Plant lifetime (years)']
        crf_battery = functions.CRF(interest_battery, lifetime_battery)
        capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} battery costs'] = \
            hexagons[f'{demand_center} {transport_method} battery capacity'] * capital_cost_battery * crf_battery
        hexagons[f'{demand_center} LCOH - {transport_method} battery costs portion'] = \
            hexagons[f'{demand_center} {transport_method} battery costs'] \
                / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # Electrolyzer
        interest_electrolyzer = country_parameters.loc[country, 'Plant interest rate']
        lifetime_electrolyzer = country_parameters.loc[country, 'Plant lifetime (years)']
        crf_electrolyzer = functions.CRF(interest_electrolyzer, lifetime_electrolyzer)
        capital_cost_electrolyzer = links_parameters.loc['Electrolysis', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} electrolyzer costs'] = \
            hexagons[f'{demand_center} {transport_method} electrolyzer capacity'] * capital_cost_electrolyzer * crf_electrolyzer
        hexagons[f'{demand_center} LCOH - {transport_method} electrolyzer portion'] = \
            hexagons[f'{demand_center} {transport_method} electrolyzer costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
 
        # H2 Storage
        interest_h2_storage = country_parameters.loc[country, 'Plant interest rate']
        lifetime_h2_storage = country_parameters.loc[country, 'Plant lifetime (years)']
        crf_h2_storage = functions.CRF(interest_h2_storage, lifetime_h2_storage)
        capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} H2 storage costs'] = \
            hexagons[f'{demand_center} {transport_method} H2 storage capacity'] * capital_cost_h2_storage * crf_h2_storage
        hexagons[f'{demand_center} LCOH - {transport_method} H2 storage portion'] = \
            hexagons[f'{demand_center} {transport_method} H2 storage costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
       
        for generator in snakemake.config['generators']:
            generator_lower = generator.lower()
            interest_generator = country_parameters.loc[country, f'{generator} interest rate']
            lifetime_generator = country_parameters.loc[country, f'{generator} lifetime (years)']
            crf_generator = functions.CRF(interest_generator, lifetime_generator)
            capital_cost_generator = generators_parameters.loc[f'{generator}', 'capital_cost']
            hexagons[f'{demand_center} {transport_method} {generator_lower} costs'] = \
                hexagons[f'{demand_center} {transport_method} {generator_lower} capacity'] * capital_cost_generator * crf_generator
            hexagons[f'{demand_center} LCOH - {transport_method} {generator_lower} portion'] = \
                hexagons[f'{demand_center} {transport_method} {generator_lower} costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

hexagons.to_file(str(snakemake.output[0]), driver='GeoJSON', encoding='utf-8')
hexagons.to_csv(str(snakemake.output[1]), encoding='latin-1')