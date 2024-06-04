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
import os
if __name__ == "__main__":
    
    logging.basicConfig(level=logging.INFO)
    
    # calculate min and max coordinates from hexagons
    hexagons = gpd.read_file(str(snakemake.input.hexagons))
    
    hexagon_bounds = hexagons.geometry.bounds
    min_lon, min_lat = hexagon_bounds[['minx','miny']].min()
    max_lon, max_lat = hexagon_bounds[['maxx','maxy']].max()

    weather_year = snakemake.wildcards.weather_year
    end_weather_year = int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    
    # Create folders for final cutouts and temporary files
    if not os.path.exists('Cutouts'):
        os.makedirs('Cutouts')
    if not os.path.exists('temp'):
        os.makedirs('temp')
    
    cutout = atlite.Cutout(
        path=str(snakemake.output),
        module="era5",
        x=slice(min_lon, max_lon),
        y=slice(min_lat, max_lat),
        time=slice(start_date, end_date),
    )
    
    cutout.prepare(tmpdir="temp") # TEMPDIR DEFINITION IS NEW TO FIX ERROR