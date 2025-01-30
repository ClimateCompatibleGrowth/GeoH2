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
from functions import CRF

# Load hexagons
hexagons = gpd.read_file(str(snakemake.input.hexagons))
generators = dict(snakemake.config['generators_dict'])

plant_type = str(snakemake.wildcards.plant_type)

# Load necessary parameters
demand_excel_path = str(snakemake.input.demand_parameters)
demand_parameters = pd.read_excel(demand_excel_path, index_col='Demand center')
country_excel_path = str(snakemake.input.country_parameters)
country_parameters = pd.read_excel(country_excel_path, index_col='Country')
if plant_type == "hydrogen":
    storage_csv_path = 'parameters/basic_h2_plant/storage_units.csv' # Battery
    storage_parameters = pd.read_csv(storage_csv_path, index_col='name')
    stores_csv_path = 'parameters/basic_h2_plant/stores.csv' # H2 storage 
    stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
    links_csv_path = 'parameters/basic_h2_plant/links.csv' # Electrolyzer 
    links_parameters = pd.read_csv(links_csv_path, index_col='name')
    generators_csv_path = 'parameters/basic_h2_plant/generators.csv' # Solar and generator 
    generators_parameters = pd.read_csv(generators_csv_path, index_col='name')
elif plant_type == "ammonia":
    stores_csv_path = 'parameters/basic_nh3_plant/stores.csv' # H2 storage 
    stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
    links_csv_path = 'parameters/basic_nh3_plant/links.csv' # Electrolyzer 
    links_parameters = pd.read_csv(links_csv_path, index_col='name')
    generators_csv_path = 'parameters/basic_nh3_plant/generators.csv' # Solar and generator 
    generators_parameters = pd.read_csv(generators_csv_path, index_col='name')


# For each demand center, get costs for each component

demand_centers = demand_parameters.index
transport_methods = ['pipeline', 'trucking']
for demand_center in demand_centers:
    print(f"\nCalculating for {demand_center} begins...")
    # Get location of demand center
    lat = demand_parameters.loc[demand_center, 'Lat [deg]']
    lon = demand_parameters.loc[demand_center, 'Lon [deg]']
    coordinates = str(lat) + ", " + str(lon)
    # Get country where the demand center is
    geolocator = Photon(user_agent="MyApp")
    location = geolocator.reverse(coordinates, language="en")
    country = location.raw['properties']['country']
    
    # Store CRF from plant data
    interest_plant = country_parameters.loc[country, 'Plant interest rate']
    lifetime_plant = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_plant = CRF(interest_plant, lifetime_plant)
    
    for transport_method in transport_methods:  
        # Work out the cost for each component using the data for the country you are looking at
        # Battery
        if plant_type == "hydrogen":
            capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
        elif plant_type == "ammonia":
            capital_cost_battery = stores_parameters.loc['Battery', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} battery costs'] = \
            hexagons[f'{demand_center} {transport_method} battery capacity'] *\
                capital_cost_battery * crf_plant
        hexagons[f'{demand_center} LC - {transport_method} battery costs portion'] = \
            hexagons[f'{demand_center} {transport_method} battery costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # Electrolyzer
        capital_cost_electrolyzer = links_parameters.loc['Electrolysis', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} electrolyzer costs'] = \
            hexagons[f'{demand_center} {transport_method} electrolyzer capacity'] *\
                capital_cost_electrolyzer * crf_plant
        hexagons[f'{demand_center} LC - {transport_method} electrolyzer portion'] = \
            hexagons[f'{demand_center} {transport_method} electrolyzer costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
 
        # H2 Storage
        if plant_type == "hydrogen":
            capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
        elif plant_type == "ammonia":
                capital_cost_h2_storage = stores_parameters.loc['CompressedH2Store', 'capital_cost']
        hexagons[f'{demand_center} {transport_method} H2 storage costs'] = \
            hexagons[f'{demand_center} {transport_method} H2 storage capacity'] *\
                capital_cost_h2_storage * crf_plant
        hexagons[f'{demand_center} LC - {transport_method} H2 storage portion'] = \
            hexagons[f'{demand_center} {transport_method} H2 storage costs']/ \
                demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
        
        # Work out CRF, then work out the cost for each generator using the data for the country you are looking at
        # Generators
        for generator in generators.keys():
            generator_lower = generator.lower()
            interest_generator = country_parameters.loc[country, f'{generator} interest rate']
            lifetime_generator = country_parameters.loc[country, f'{generator} lifetime (years)']
            crf_generator = CRF(interest_generator, lifetime_generator)
            capital_cost_generator = generators_parameters.loc[f'{generator}', 'capital_cost']
            hexagons[f'{demand_center} {transport_method} {generator_lower} costs'] = \
                hexagons[f'{demand_center} {transport_method} {generator_lower} capacity'] * capital_cost_generator * crf_generator
            hexagons[f'{demand_center} LC - {transport_method} {generator_lower} portion'] = \
                hexagons[f'{demand_center} {transport_method} {generator_lower} costs']/ \
                    demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

print("\nCalculations complete.\n")
hexagons.to_file(snakemake.output[0], driver='GeoJSON', encoding='utf-8') # snakemake config
hexagons.to_csv(snakemake.output[1], encoding='latin-1') # snakemake config