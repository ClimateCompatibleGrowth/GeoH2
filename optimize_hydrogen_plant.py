# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:47:57 2023

@author: Claire Halloran, University of Oxford

Includes code from Nicholas Salmon, University of Oxford, for optimizing
hydrogen plant capacity.

"""

import atlite
import geopandas as gpd
import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import p_H2_aux as aux
from functions import NPV
import numpy as np
import logging
logging.basicConfig(level=logging.ERROR)

hexagons = gpd.read_file('Data/hex_total.geojson')

#!!! currently ignoring land-use restrictions-- add later
# excluder = atlite.gis.ExclusionContainer()

cutout = atlite.Cutout("Cutouts/Kenya-2022-07.nc")

layout = cutout.uniform_layout()
# can add hydro layout here if desired using hydrogen potential map

pv_profile = cutout.pv(
    panel= 'CSi',
    orientation='latitude_optimal',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
pv_profile = pv_profile.rename(dict(dim_0='hexagon'))

wind_profile = cutout.wind(
    turbine = 'Vestas_V80_2MW_gridstreamer',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
wind_profile = wind_profile.rename(dict(dim_0='hexagon'))

def optimize_hydrogen_plant(wind_potential, pv_potential, times, interest, basis_fn = None):
    # ==================================================================================================================
    # Set up network
    # ==================================================================================================================
    # print(wind_potential.shape)
    # print(pv_potential.shape)
    # Import a generic network
    n = pypsa.Network(override_component_attrs=aux.create_override_components())

    # Set the time values for the network
    # n.set_snapshots(range(len(weather_data)))
    n.set_snapshots(times)

    # Import the design of the H2 plant into the network
    n.import_from_csv_folder("Data/Basic_H2_plant")

    # renewables_list = n.generators.index.to_list()
    # Note: All flows are in MW or MWh, conversions for hydrogen done using HHVs. Hydrogen HHV = 39.4 MWh/t

    # ==================================================================================================================
    # Send the weather data to the model
    # ==================================================================================================================
    #!!! need to implement a max capacity based on land availability
    n.generators_t.p_max_pu['Wind'] = wind_potential
    n.generators_t.p_max_pu['Solar'] = pv_potential

    # ==================================================================================================================
    # Check if the CAPEX input format in Basic_H2_plant is correct, and fix it up if not
    # ==================================================================================================================
    lifetime = 20 #!!! temporary assumption of generation asset lifetime, add to inputs
    # CAPEX_check = aux.check_CAPEX()
    # !!! need to specify technology-specific and country-specific WACC here
    for item in [n.generators, n.links, n.stores]:
        item.capital_cost = item.capital_cost * NPV(interest,lifetime)

    # ==================================================================================================================
    # Solve the model
    # ==================================================================================================================

    # Ask the user how they would like to solve
    solver = 'cbc'
    if basis_fn is None:
        n.lopf(solver_name=solver,
               pyomo=False,
               extra_functionality=aux.extra_functionalities,
               store_basis = True
               )
    else:
        n.lopf(solver_name=solver,
               pyomo=False,
               extra_functionality=aux.extra_functionalities,
               warmstart = basis_fn,
               store_basis = True
               )
    # ==================================================================================================================
    # Output results
    # ==================================================================================================================
    # output = {
    #         #!!! need to redo this with a load profile
    #     'LCOH (USD/kg)': n.objective/(n.loads.p_set.values[0]/39.4*8760*1000),
    #     'Wind capacity (MW)': n.generators.p_nom_opt['Wind'],
    #     'Solar capacity (MW)': n.generators.p_nom_opt['Solar'],
    #     'Electrolyzer capacity (MW)': n.links.p_nom_opt['Electrolysis'],
    #     'H2 storage capacity (MWh)': n.stores.e_nom_opt['Compressed H2 Store'],
    # }
    # output_array = np.array([
    #         #!!! need to redo this with a load profile
    #     n.objective/(n.loads.p_set.values[0]/39.4*8760*1000),
    #     n.generators.p_nom_opt['Wind'],
    #     n.generators.p_nom_opt['Solar'],
    #     n.links.p_nom_opt['Electrolysis'],
    #     n.stores.e_nom_opt['Compressed H2 Store'],
    # ])
    #!!! for now just return LCOH
    lcoh = n.objective/(n.loads.p_set.values[0]/39.4*8760*1000)
    print(lcoh)
    basis_fn = n.basis_fn
    return lcoh, basis_fn # need to change this to something that can be added to hexagons



# test run-- run a single hexagon to get warmstart
# test = optimize_hydrogen_plant(wind_profile.sel(hexagon=1),pv_profile.sel(hexagon=1),pv_profile.time, 0.08)

# test = optimize_hydrogen_plant(wind_profile.sel(hexagon=hexagon_number),
#                             pv_profile.sel(hexagon=hexagon_number),
#                             0.08)
# info to save
# production cost in USD or euros per kg
# wind and solar optimal capacity
# electrolyzer size


# chunk hexagons for dask parallelization

# wind_profile = wind_profile.chunk({'hexagon':1})
# pv_profile = pv_profile.chunk({'hexagon':1})

#%% try speeding up with a warmstart for loop

import time

interest_rate = 0.08
lcohs = np.zeros(len(pv_profile.hexagon))
bases = None

start = time.process_time()

for hexagon in pv_profile.hexagon.data:
    if bases == None:
        lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                pv_profile.sel(hexagon = hexagon),
                                wind_profile.time,
                                interest_rate, 
                                # !!! add a hexagon-based interest rate
                                )
    else:
        print('Warmstarting...')
        lcoh, basis_fn = optimize_hydrogen_plant(wind_profile.sel(hexagon = hexagon),
                                pv_profile.sel(hexagon = hexagon),
                                wind_profile.time,
                                interest_rate,
                                # !!! add a hexagon-based interest rate
                                basis_fn = bases
                                )
    lcohs[hexagon]=lcoh
    bases = basis_fn
no_compute = time.process_time()-start
print(str(no_compute) + ' s')

#%% add optimal LCOH for each hexagon to hexagon file
hexagons['LCOH'] = lcohs

hexagons.to_file('Resources/hex_lcoh.geojson')
#%% plot LCOH for each hexagon
# update central coordinates for area considered
crs = ccrs.Orthographic(central_longitude = 37.5, central_latitude= 0.0)

fig = plt.figure(figsize=(10,5))

ax = plt.axes(projection=crs)
ax.set_axis_off()

hexagons.to_crs(crs.proj4_init).plot(
    ax=ax,
    column = 'LCOH',
    legend = True,
    cmap = 'viridis_r',
    legend_kwds={'label':'Production LCOH [$/kg]'}
    
)
