# -*- coding: utf-8 -*-
"""
This module...


"""
import geopandas as gpd
import modify_hex_file as modify
import get_weather_data as weather
import json


def main():
    # inputs and arguments
    hexagon_path = "data/hex_final.geojson" # SNAKEMAKE INPUT
    country_parameters = "parameters/country_parameters.xlsx" # SNAKEMAKE INPUT
    hexagons = gpd.read_file(hexagon_path) 
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')) # may need to switch to higher res
    
    # updating hexagon file
    hexagons_with_country = modify.assign_country(hexagons, world)
    hexagons_with_country.to_file("data/updated_hex_final.geojson", driver="GeoJSON")
    updated_hexagons = modify.drop_extra_hexagons(hexagon_path, country_parameters)
    with open(hexagon_path, 'w') as file:
        json.dump(updated_hexagons, file)
   
    # storing weather data
    weather.get_weather_date(hexagons)

if __name__ == "__main__":
    main()