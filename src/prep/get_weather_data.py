# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 15:11:47 2023

@author: Claire Halloran, University of Oxford

get_weather_data.py

This script fetches historical weather data to calculate wind and solar potential
from the ERA-5 reanalysis dataset using Atlite.

Create cutouts with `atlite <https://atlite.readthedocs.io/en/latest/>`_.

For this rule to work you must have

- installed the `Copernicus Climate Data Store <https://cds.climate.copernicus.eu>`_ ``cdsapi`` package  (`install with `pip``) and
- registered and setup your CDS API key as described on their website <https://cds.climate.copernicus.eu/api-how-to>

"""
import logging
import atlite
import geopandas as gpd
from utils import check_folder_exists

def calculate_coords(hexagons):
    """
    Calculates mininum and maximum coordinates using bounds from the hexagon file.

    ...
    Parameters
    ----------
    hexagons : geodataframe
        Updated hexagon file from data folder.

    Returns
    -------
    min_lon : int
        Mininum longitude
    min_lat : int
        Minimum latitude
    max_lon : int
        Maximum longitude
    max_lat : int
        Maximum latitude
    """
    hexagon_bounds = hexagons.geometry.bounds
    min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
    max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()
    
    return min_lon, min_lat, max_lon, max_lat

def prepare_cutout(min_lon, min_lat, max_lon, max_lat, start_date, end_date):
    """
    Creates and prepares the cutout, using Atlite.

    ...
    Parameters
    ----------
    min_lon : int
        Mininum longitude
    min_lat : int
        Minimum latitude
    max_lon : int
        Maximum longitude
    max_lat : int
        Maximum latitude
    start_date : string
        Start date for weather collection in 'YYYY-MM-DD' format
    end_date : string
        End date for weather collection in 'YYYY-MM-DD' format
    """
    cutout = atlite.Cutout(
        path=str(snakemake.output),
        module="era5",
        x=slice(min_lon, max_lon),
        y=slice(min_lat, max_lat),
        time=slice(start_date, end_date),
    )
    
    cutout.prepare(tmpdir="temp", show_progress=True) # TEMPDIR DEFINITION IS NEW TO FIX ERROR

def main():
    try:
        hexagons = gpd.read_file(f"data/hexagons_with_country_{snakemake.wildcards.country}_{snakemake.config['scenario']['plant_type']}.geojson")
    except:
        print(f"There is no file called \'hexagons_with_country_{snakemake.wildcards.country}_{snakemake.config['scenario']['plant_type']}.geojson\' in the 'data/' folder. \nRun the necessary rule to produce that file.")
    else:
        # Displays information on process as it runs.
        logging.basicConfig(level=logging.INFO)

        min_lon, min_lat, max_lon, max_lat = calculate_coords(hexagons)

        start_weather_year = int(snakemake.wildcards.weather_year)
        end_weather_year = int(snakemake.wildcards.weather_year)+int(snakemake.config["years_to_check"])
        start_date = f'{start_weather_year}-01-01'
        end_date = f'{end_weather_year}-01-01'

        check_folder_exists("cutouts")
        check_folder_exists("temp")

        prepare_cutout(min_lon, min_lat, max_lon, max_lat, start_date, end_date)

if __name__ == "__main__":
    main()