"""
Created on Tue Aug 29th 2023

@author: Alycia Leonard, University of Oxford, alycia.leonard@eng.ox.ac.uk

Cost by hydrogen plant component

Add attributes to hex file for cost of each component
"""

# %% identify lowest-cost strategy: trucking vs. pipeline

import geopandas as gpd
import pandas as pd
import numpy as np
from geopy.geocoders import Photon
import functions

# Load hexagons
hexagons = gpd.read_file('Resources/hex_total_cost.geojson')

# Load necessary parameters
demand_excel_path = 'Parameters/demand_parameters.xlsx'
demand_parameters = pd.read_excel(demand_excel_path, index_col='Demand center')
country_excel_path = 'Parameters/country_parameters.xlsx'
country_parameters = pd.read_excel(country_excel_path, index_col='Country')
stores_csv_path = 'Parameters/Basic_H2_plant/stores.csv' # H2 storage
stores_parameters = pd.read_csv(stores_csv_path, index_col='name')
storage_csv_path = 'Parameters/Basic_H2_plant/storage_units.csv' # Battery
storage_parameters = pd.read_csv(storage_csv_path, index_col='name')
links_csv_path = 'Parameters/Basic_H2_plant/links.csv' # Electrolyzer
links_parameters = pd.read_csv(links_csv_path, index_col='name')
generators_csv_path = 'Parameters/Basic_H2_plant/generators.csv' # Solar and wind
generators_parameters = pd.read_csv(generators_csv_path, index_col='name')

# For each demand center, get costs for each component

demand_centers = demand_parameters.index
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

    # Battery - pipeline 
    interest_battery = country_parameters.loc[country, 'Plant interest rate']
    lifetime_battery = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_battery = functions.CRF(interest_battery, lifetime_battery)
    capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
    hexagons[f'{demand_center} pipeline battery costs'] = \
        hexagons[f'{demand_center} pipeline battery capacity'] * capital_cost_battery * crf_battery
    hexagons[f'{demand_center} LCOH - pipeline battery costs portion'] = \
        hexagons[f'{demand_center} pipeline battery costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
    
    # Battery - trucking 
    interest_battery = country_parameters.loc[country, 'Plant interest rate']
    lifetime_battery = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_battery = functions.CRF(interest_battery, lifetime_battery)
    capital_cost_battery = storage_parameters.loc['Battery', 'capital_cost']
    hexagons[f'{demand_center} trucking battery costs'] = \
        hexagons[f'{demand_center} trucking battery capacity'] * capital_cost_battery * crf_battery
    hexagons[f'{demand_center} LCOH - trucking battery costs portion'] = \
        hexagons[f'{demand_center} trucking battery costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # Electrolyzer - pipeline
    interest_electrolyzer = country_parameters.loc[country, 'Plant interest rate']
    lifetime_electrolyzer = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_electrolyzer = functions.CRF(interest_electrolyzer, lifetime_electrolyzer)
    capital_cost_electrolyzer = links_parameters.loc['Electrolysis', 'capital_cost']
    hexagons[f'{demand_center} pipeline electrolyzer costs'] = \
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] * capital_cost_electrolyzer * crf_electrolyzer
    hexagons[f'{demand_center} LCOH - pipeline electrolyzer portion'] = \
        hexagons[f'{demand_center} pipeline electrolyzer costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # Electrolyzer - trucking
    interest_electrolyzer = country_parameters.loc[country, 'Plant interest rate']
    lifetime_electrolyzer = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_electrolyzer = functions.CRF(interest_electrolyzer, lifetime_electrolyzer)
    capital_cost_electrolyzer = links_parameters.loc['Electrolysis', 'capital_cost']
    hexagons[f'{demand_center} trucking electrolyzer costs'] = \
        hexagons[f'{demand_center} trucking electrolyzer capacity'] * capital_cost_electrolyzer * crf_electrolyzer
    hexagons[f'{demand_center} LCOH - trucking electrolyzer portion'] = \
        hexagons[f'{demand_center} trucking electrolyzer costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # H2 Storage - pipeline
    interest_h2_storage = country_parameters.loc[country, 'Plant interest rate']
    lifetime_h2_storage = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_h2_storage = functions.CRF(interest_h2_storage, lifetime_h2_storage)
    capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
    hexagons[f'{demand_center} pipeline H2 storage costs'] = \
        hexagons[f'{demand_center} pipeline H2 storage capacity'] * capital_cost_h2_storage * crf_h2_storage
    hexagons[f'{demand_center} LCOH - pipeline H2 storage portion'] = \
        hexagons[f'{demand_center} pipeline H2 storage costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # H2 Storage - trucking
    interest_h2_storage = country_parameters.loc[country, 'Plant interest rate']
    lifetime_h2_storage = country_parameters.loc[country, 'Plant lifetime (years)']
    crf_h2_storage = functions.CRF(interest_h2_storage, lifetime_h2_storage)
    capital_cost_h2_storage = stores_parameters.loc['Compressed H2 Store', 'capital_cost']
    hexagons[f'{demand_center} trucking H2 storage costs'] = \
        hexagons[f'{demand_center} trucking H2 storage capacity'] * capital_cost_h2_storage * crf_h2_storage
    hexagons[f'{demand_center} LCOH - trucking H2 storage portion'] = \
        hexagons[f'{demand_center} trucking H2 storage costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # Wind - pipeline
    interest_wind = country_parameters.loc[country, 'Wind interest rate']
    lifetime_wind = country_parameters.loc[country, 'Wind lifetime (years)']
    crf_wind = functions.CRF(interest_wind, lifetime_wind)
    capital_cost_wind = generators_parameters.loc['Wind', 'capital_cost']
    hexagons[f'{demand_center} pipeline wind costs'] = \
        hexagons[f'{demand_center} pipeline wind capacity'] * capital_cost_wind * crf_wind
    hexagons[f'{demand_center} LCOH - pipeline wind portion'] = \
        hexagons[f'{demand_center} pipeline wind costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
    
    # Wind - trucking
    interest_wind = country_parameters.loc[country, 'Wind interest rate']
    lifetime_wind = country_parameters.loc[country, 'Wind lifetime (years)']
    crf_wind = functions.CRF(interest_wind, lifetime_wind)
    capital_cost_wind = generators_parameters.loc['Wind', 'capital_cost']
    hexagons[f'{demand_center} trucking wind costs'] = \
        hexagons[f'{demand_center} trucking wind capacity'] * capital_cost_wind * crf_wind
    hexagons[f'{demand_center} LCOH - trucking wind portion'] = \
        hexagons[f'{demand_center} trucking wind costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

    # Solar - pipeline
    interest_solar = country_parameters.loc[country, 'Solar interest rate']
    lifetime_solar = country_parameters.loc[country, 'Solar lifetime (years)']
    crf_solar = functions.CRF(interest_solar, lifetime_solar)
    capital_cost_solar = generators_parameters.loc['Solar', 'capital_cost']
    hexagons[f'{demand_center} pipeline solar costs'] = \
        hexagons[f'{demand_center} pipeline solar capacity'] * capital_cost_solar * crf_solar
    hexagons[f'{demand_center} LCOH - pipeline solar portion'] = \
        hexagons[f'{demand_center} pipeline solar costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']
    
    # Solar - trucking
    interest_solar = country_parameters.loc[country, 'Solar interest rate']
    lifetime_solar = country_parameters.loc[country, 'Solar lifetime (years)']
    crf_solar = functions.CRF(interest_solar, lifetime_solar)
    capital_cost_solar = generators_parameters.loc['Solar', 'capital_cost']
    hexagons[f'{demand_center} trucking solar costs'] = \
        hexagons[f'{demand_center} trucking solar capacity'] * capital_cost_solar * crf_solar
    hexagons[f'{demand_center} LCOH - trucking solar portion'] = \
        hexagons[f'{demand_center} trucking solar costs'] / demand_parameters.loc[demand_center, 'Annual demand [kg/a]']

hexagons.to_file('Resources/hex_cost_components.geojson', driver='GeoJSON', encoding='utf-8')
hexagons.to_csv('Resources/hex_cost_components.csv', encoding='latin-1')