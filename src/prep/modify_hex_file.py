# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 12:05:51 2023

@author: Claire Halloran, University of Oxford

Functions for preparing the hexagon file for optimization.

"""
import geopandas as gpd
import pandas as pd
import json
import warnings
warnings.filterwarnings("ignore")

def assign_country(hexagons, world):
    """
    Assigns interest rate to different hexagons for different 
    technology categories based on their country.

    ...
    Parameters
    ----------
    hexagons : geodataframe
        Hexagon file from data folder
    world : geodataframe
        World dataset

    Returns
    -------
    hexagons_with_country : geodataframe
        Modified hexagons
    """
    hexagons.to_crs(world.crs, inplace=True)
    countries = world.drop(columns=[
                                    'pop_est', 
                                    'continent', 
                                    'iso_a3', 
                                    'gdp_md_est',
                                    ]
                            )
    countries = countries.rename(columns={'name':'country'})
    hexagons_with_country = gpd.sjoin(hexagons, countries, op='intersects') # changed from "within"

    return hexagons_with_country

def remove_extra_hexagons(hexagon_path, country_parameters):
    """
    Removes duplicated hexagons.

    ...
    Parameters
    ----------
    hexagon_path : string
        File path to hexagon file
    country_parameters : dataframe
        Country file from parameters

    Returns
    -------
    hexagons : geodataframe
        Modified hexagons
    """
    with open(hexagon_path, 'r') as file:
        hexagons = json.load(file)

    copied_list = hexagons["features"].copy()
    country_name = country_parameters.index.values[0]

    for feature in copied_list:
        if feature['properties']['country'] != country_name:
            hexagons['features'].remove(feature)

    return hexagons

def update_hexagons(hexagons, hexagon_path):
    """
    Updates hexagon file with the new information
    """
    if isinstance(hexagons, dict):
        with open(hexagon_path, 'w') as file:
            json.dump(hexagons, file)
    else:
        hexagons.to_file(f"{hexagon_path}", driver="GeoJSON")

def main():
    hexagon_path = "data/hex_final_DJ.geojson" # SNAKEMAKE INPUT
    country_parameters = pd.read_excel("parameters/country_parameters.xlsx",
                                        index_col='Country') # SNAKEMAKE INPUT
    hexagons = gpd.read_file(hexagon_path) 
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')) # may need to switch to higher res

    hexagons_with_country = assign_country(hexagons, world)
    update_hexagons(hexagons_with_country, hexagon_path)

    # Finish off with the removing extra hexagons.
    final_hexagons = remove_extra_hexagons(hexagon_path, country_parameters)
    update_hexagons(final_hexagons, hexagon_path)


if __name__ == "__main__":
    main()