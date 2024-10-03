# -*- coding: utf-8 -*-
"""
This module...


"""
import geopandas as gpd
import json


def assign_country(hexagons, world):
    """
    Assigns interest rate to different hexagons for different 
    technology categories based on their country.

    ...
    Parameters
    ----------
    hexagons : geodataframe
        Hexagon file from data folder.
    world : geodataframe
        World dataset.

    Returns
    -------
    hexagons_with_country : geodataframe
        Modified hexagons.
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

def drop_extra_hexagons(hexagon_path, country_parameters):
    """
    Removes duplicated hexagons.

    ...
    Parameters
    ----------
    hexagons : geodataframe
        Geodataframe containing the hex file.
    """
    with open(hexagon_path, 'r') as file:
        hexagons = json.load(file)

    copied_list = hexagons["features"].copy()

    for feature in copied_list:
        if feature['properties']['country'] != country_parameters.index.values[0]:
            hexagons['features'].remove(feature)

    return hexagons

def update_hexagons(hexagons, hexagon_path=None):
    """
    Updates hexagon file with the new information
    """
    if isinstance(hexagons, dict):
        with open(hexagon_path, 'w') as file:
            json.dump(hexagons, file)
    else:
        hexagons.to_file("data/updated_hex_final.geojson", driver="GeoJSON")