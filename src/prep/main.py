# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 12:05:51 2023

@author: Claire Halloran, University of Oxford

Functions for preparing the hexagon file for optimization.

"""
import geopandas as gpd
import pandas as pd
import json

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

def remove_extra_hexagons(output_hexagon_path, country_parameters):
    """
    Removes duplicated hexagons.

    ...
    Parameters
    ----------
    output_hexagon_path : string
        File path to output hexagon file
    country_parameters : dataframe
        Country file from parameters

    Returns
    -------
    hexagons : geodataframe
        Modified hexagons
    """
    with open(output_hexagon_path, 'r') as file:
        hexagons = json.load(file)

    copied_list = hexagons["features"].copy()
    country_name = country_parameters.index.values[0]

    for feature in copied_list:
        if feature['properties']['country'] != country_name:
            hexagons['features'].remove(feature)

    return hexagons

def update_hexagons(hexagons, output_hexagon_path):
    """
    Updates hexagon file with the new information
    """
    if isinstance(hexagons, dict):
        with open(output_hexagon_path, 'w') as file:
            json.dump(hexagons, file)
    else:
        hexagons.to_file(f"{output_hexagon_path}", driver="GeoJSON")

def main():
    print("Prepping file...")
    country_parameters = pd.read_excel(str(snakemake.input.country_parameters),
                                        index_col='Country')
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')) # may need to switch to higher res

    output_hexagon_path = str(snakemake.output)

    hexagons_with_country = assign_country(hexagons, world)
    update_hexagons(hexagons_with_country, output_hexagon_path)

    # Finish off with the removing extra hexagons.
    final_hexagons = remove_extra_hexagons(output_hexagon_path, country_parameters)
    update_hexagons(final_hexagons, output_hexagon_path)
    print("File prepped.")

if __name__ == "__main__":
    main()