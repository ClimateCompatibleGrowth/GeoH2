#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 12:05:51 2023

@author: Claire Halloran, University of Oxford

Assigns interest rate to different hexagons for different technology categories
based on their country.

Just add to optimize_hydrogen_plant.py?

"""
import geopandas as gpd

hexagons = gpd.read_file('Data/hex_final.geojson')
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')) # may need to switch to higher res
countries = world.drop(columns=['pop_est', 'continent', 'iso_a3', 'gdp_md_est'])
countries = countries.rename(columns={'name':'country'})
hexagons_with_country = gpd.sjoin(hexagons, countries, op='within')
hexagons_with_country.to_file('Data/hexagons_with_country.geojson', driver='GeoJSON')