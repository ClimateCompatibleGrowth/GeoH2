# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:47:57 2023

@author: Claire Halloran, University of Oxford

Includes code from Nicholas Salmon, University of Oxford, for optimizing
hydrogen plant capacity.

"""

import atlite
from functions import CRF
import geopandas as gpd
import logging
import numpy as np
import pandas as pd
import pypsa
import warnings


class Plant:
    def __init__(self, n, type, water_limit=None):
        self.n = n
        self.type = type
        self.water_limit = water_limit

    def get_hydrogen_network(self, demand_profile, times, country_series):
        '''
        Sets up the network.

        Parameters
        ----------
        demand_profile : pandas DataFrame
            hourly dataframe of hydrogen demand in kg.
        times : xarray DataArray
            1D dataarray with timestamps for wind and solar potential.
        country_series: pandas Series
            interest rate and lifetime information.

        Returns
        -------
        n :
            pypsa network.
        '''
        # Set the time values for the network
        self.n.set_snapshots(times)

        # Import the design of the H2 plant into the network
        self.n.import_from_csv_folder("parameters/basic_h2_plant")

        # Import demand profile
        # Note: All flows are in MW or MWh, conversions for hydrogen done using HHVs. Hydrogen HHV = 39.4 MWh/t
        self.n.add('Load',
            'Hydrogen demand',
            bus = 'Hydrogen',
            p_set = demand_profile['Demand']/1000*39.4,
            )

        # need to discuss with Alycia, does this need it's own function, would it be the same for any plant? 
        # should it be in the main body?
        for item in [self.n.links, self.n.stores, self.n.storage_units]:
            item.capital_cost = item.capital_cost * CRF(country_series['Plant interest rate'],country_series['Plant lifetime (years)'])
        
        return n
    
    # think about the below function and if it makes sense
    def set_generator_in_network(self, generator, generator_potential, generator_max_capacity, country_series):
        '''
        Sets provided generator in the network.

        Parameters
        ----------
        generator : string
            name of generator to add into the network.
        generator_potential : xarray DataArray
            1D dataarray of per-unit generator potential in hexagon.
        generator_max_capacity :
            ...
        country_series: pandas Series
            interest rate and lifetime information.
        '''
        # Send the weather data to the model
        self.n.generators_t.p_max_pu[f'{generator}'] = generator_potential

        # specify maximum capacity based on land use
        self.n.generators.loc[f'{generator}','p_nom_max'] = generator_max_capacity*4

        # specify technology-specific and country-specific WACC and lifetime here
        self.n.generators.loc[f'{generator}','capital_cost'] = self.n.generators.loc[f'{generator}','capital_cost']\
            * CRF(country_series[f'{generator} interest rate'], country_series[f'{generator} lifetime (years)'])
        
    def set_wind_in_network(self, wind_potential, wind_max_capacity, country_series):
        '''
        Sets wind in the network.

        Parameters
        ----------
        generator : string
            name of generator to add into the network.
        wind_potential : xarray DataArray
            1D dataarray of per-unit wind potential in hexagon.
        wind_max_capacity :
            ...
        country_series: pandas Series
            interest rate and lifetime information.
        '''
         # Send the weather data to the model
        self.n.generators_t.p_max_pu['Wind'] = wind_potential

        # specify maximum capacity based on land use
        self.n.generators.loc['Wind','p_nom_max'] = wind_max_capacity*4

        # specify technology-specific and country-specific WACC and lifetime here
        self.n.generators.loc['Wind','capital_cost'] = self.n.generators.loc['Wind','capital_cost']\
            * CRF(country_series['Wind interest rate'], country_series['Wind lifetime (years)'])
        
    def set_pv_in_network(self, pv_potential, pv_max_capacity, country_series):
        '''
        Sets solar in the network.

        Parameters
        ----------
        generator : string
            name of generator to add into the network.
        pv_potential : xarray DataArray
            1D dataarray of per-unit solar potential in hexagon.
        pv_max_capacity :
            ...
        country_series: pandas Series
            interest rate and lifetime information.
        '''
         # Send the weather data to the model
        self.n.generators_t.p_max_pu['Solar'] = pv_potential

        # specify maximum capacity based on land use
        self.n.generators.loc['Solar','p_nom_max'] = pv_max_capacity

        # specify technology-specific and country-specific WACC and lifetime here
        self.n.generators.loc['Solar','capital_cost'] = self.n.generators.loc['Solar','capital_cost']\
            * CRF(country_series['Solar interest rate'], country_series['Solar lifetime (years)'])

    def solve_model(self, solver):
        '''
        Solves model using the provided solver.

        Parameters
        ----------
        solver : string
            name of solver to be used.
        '''
        self.n.lopf(solver_name=solver,
            solver_options = {'LogToConsole':0, 'OutputFlag':0},
            pyomo=False,
            )

    def get_water_constraint(self, demand_profile): 
        '''
        Calculates the water constraint.

        Parameters
        ----------
        demand_profile : pandas DataFrame
            hourly dataframe of hydrogen demand in kg.
        water_limit : float
            annual limit on water available for electrolysis in hexagon, in cubic meters.

        Returns
        -------
        water_constraint : boolean
            whether there is a water constraint or not.
        '''
        # total hydrogen demand in kg
        total_hydrogen_demand = demand_profile['Demand'].sum()
        # check if hydrogen demand can be met based on hexagon water availability
        water_constraint =  total_hydrogen_demand <= self.water_limit * 111.57 # kg H2 per cubic meter of water
        
        return water_constraint

    def get_results(self, generators):
        '''
        Calculates the water constraint.

        Parameters
        ----------
        generators : list
            contains types of generators that this plant uses.
        water_limit : float
            annual limit on water available for electrolysis in hexagon, in cubic meters. Default is None.

        Returns
        -------
        lcoh : float
            levelized cost per kg hydrogen.
        generator_capactities : dictionary
            contains each generator with their optimal capacity in MW.
        electrolyzer_capacity: float
            optimal electrolyzer capacity in MW.
        battery_capacity: float
            optimal battery storage capacity in MW/MWh (1 hour batteries).
        h2_storage: float
            optimal hydrogen storage capacity in MWh.
        '''
        generator_capacities = {}
        if self.water_limit != None :
            water_constraint = plant.get_water_constraint(hydrogen_demand_trucking)
            if water_constraint == False:
                        print('Not enough water to meet hydrogen demand!')
                        lcoh = np.nan
                        for generator in generators:
                             generator_capacities[generator] = np.nan
                        electrolyzer_capacity = np.nan
                        battery_capacity = np.nan
                        h2_storage = np.nan

        if self.water_limit == None :
            lcoh = self.n.objective/(self.n.loads_t.p_set.sum()[0]/39.4*1000) # convert back to kg H2
            for generator in generators:
                 generator_capacities[generator] = self.n.generators.p_nom_opt[f"{generator}"]
            electrolyzer_capacity = self.n.links.p_nom_opt['Electrolysis']
            battery_capacity = self.n.storage_units.p_nom_opt['Battery']
            h2_storage = self.n.stores.e_nom_opt['Compressed H2 Store']

        return lcoh, generator_capacities, electrolyzer_capacity, battery_capacity, h2_storage

def create_override_components():
    """Set up new component attributes as required"""
    # Modify the capacity of a link so that it can attach to 2 buses.
    override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )

    override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
    override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]
    return override_component_attrs

def get_pv_profile(cutout, layout, hexagons):
    '''
    Sets the solar profile in the cutout.

    Parameters
    ----------
    cutout : 
        ...
    layout : 
        ...
    hexagons :
        ...
    Returns
    -------
    pv_profile : 
        ...
    '''
    pv_profile = cutout.pv(
    panel= 'CSi',
    orientation='latitude_optimal',
    layout = layout,
    shapes = hexagons,
    per_unit = True
    )
    pv_profile = pv_profile.rename(dict(dim_0='hexagon'))
    return pv_profile

def get_wind_profile(cutout, layout, hexagons):
    '''
    Sets the wind profile in the cutout.

    Parameters
    ----------
    cutout : 
        ...
    layout : 
        ...
    hexagons :
        ...
    Returns
    -------
    wind_profile : 
        ...
    '''
    wind_profile = cutout.wind(
        turbine = 'NREL_ReferenceTurbine_2020ATB_4MW',
        layout = layout,
        shapes = hexagons,
        per_unit = True
        )
    wind_profile = wind_profile.rename(dict(dim_0='hexagon'))
    return wind_profile

def get_demand_schedule(quantity, start_date, end_date, transport_state, transport_params_filepath):
    '''
    Calculates hourly hydrogen demand for truck shipment and pipeline transport.

    Parameters
    ----------
    quantity : float
        annual amount of hydrogen to transport in kilograms.
    start_date: string
        start date for demand schedule in the format YYYY-MM-DD.
    end_date: string
        end date for demand schedule in the format YYYY-MM-DD.
    transport_state : string
        state hydrogen is transported in, one of '500 bar', 'LH2', 'LOHC', or 'NH3'.
    transport_params_filepath : string
        path to transport_parameters.xlsx file

    Returns
    -------
    trucking_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for hydrogen trucking.
    pipeline_hourly_demand_schedule : pandas DataFrame
        hourly demand profile for pipeline transport.
    '''
    # schedule for pipeline
    index = pd.date_range(start_date, end_date, freq = 'H')
    pipeline_hourly_quantity = quantity/index.size
    pipeline_hourly_demand_schedule = pd.DataFrame(pipeline_hourly_quantity, index=index,  columns = ['Demand'])

    # if demand center is in hexagon
    if transport_state=="None":
        # schedule for trucking
        annual_deliveries = 365*24
        trucking_hourly_demand = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(trucking_hourly_demand, index=index, columns = ['Demand'])
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)
    else:
        transport_params = pd.read_excel(transport_params_filepath,
                                            sheet_name = transport_state,
                                            index_col = 'Parameter'
                                            ).squeeze('columns')

        truck_capacity = transport_params['Net capacity (kg H2)']

        # schedule for trucking
        annual_deliveries = quantity/truck_capacity
        quantity_per_delivery = quantity/annual_deliveries
        index = pd.date_range(start_date, end_date, periods=annual_deliveries)
        trucking_demand_schedule = pd.DataFrame(quantity_per_delivery, index=index, columns=['Demand'])
        trucking_hourly_demand_schedule = trucking_demand_schedule.resample('H').sum().fillna(0.)

    return trucking_hourly_demand_schedule, pipeline_hourly_demand_schedule


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    logging.basicConfig(level=logging.ERROR)

    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')
    country_series = country_params.iloc[0]
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center',
                                    )
    demand_centers = demand_center_list.index

    weather_year = 2022             # snakemake.wildcards.weather_year
    end_weather_year = 2023         # int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = "gurobi" # maybe make this into a snakemake wildcard?
    generators = ["Wind", "Solar"] # already in the config, used in map_costs.py line 268
    hexagons = gpd.read_file('results/completed_hex_DJ.geojson') # SNAKEMAKE INPUT

    cutout = atlite.Cutout('cutouts/DJ_2022.nc') # SNAKEMAKE INPUT
    layout = cutout.uniform_layout()
    
    pv_profile = get_pv_profile(cutout, layout, hexagons)
    wind_profile = get_wind_profile(cutout, layout, hexagons)

    for demand_center in demand_centers:
        # trucking variables
        lcohs_trucking = np.zeros(len(pv_profile.hexagon))
        t_solar_capacities= np.zeros(len(pv_profile.hexagon))
        t_wind_capacities= np.zeros(len(pv_profile.hexagon))
        t_electrolyzer_capacities= np.zeros(len(pv_profile.hexagon))
        t_battery_capacities = np.zeros(len(pv_profile.hexagon))
        t_h2_storages= np.zeros(len(pv_profile.hexagon))
        
        # pipeline variables
        lcohs_pipeline = np.zeros(len(pv_profile.hexagon))
        p_solar_capacities= np.zeros(len(pv_profile.hexagon))
        p_wind_capacities= np.zeros(len(pv_profile.hexagon))
        p_electrolyzer_capacities= np.zeros(len(pv_profile.hexagon))
        p_battery_capacities = np.zeros(len(pv_profile.hexagon))
        p_h2_storages= np.zeros(len(pv_profile.hexagon))

        hydrogen_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']

        for i in range(len(hexagons)):
            trucking_state = hexagons.loc[i, f'{demand_center} trucking state']
            wind_potential = wind_profile.sel(hexagon = i)
            pv_potential = pv_profile.sel(hexagon = i)
            wind_max_capacity = hexagons.loc[i,'theo_turbines']
            pv_max_capacity = hexagons.loc[i,'theo_pv']

            hydrogen_demand_trucking, hydrogen_demand_pipeline =\
                get_demand_schedule(hydrogen_quantity,
                                start_date,
                                end_date,
                                trucking_state,
                                transport_params_filepath)
            
            transport_types = ["trucking", "pipeline"]
            for transport in transport_types:
                n = pypsa.Network(override_component_attrs=create_override_components())
                plant = Plant(n, "Hydrogen")
                if transport == "trucking":
                    n = plant.get_hydrogen_network(hydrogen_demand_trucking, wind_profile.time, country_series)
                    plant.set_wind_in_network(wind_potential, wind_max_capacity, country_series)
                    plant.set_pv_in_network(pv_potential, pv_max_capacity, country_series)
                    plant.solve_model(solver)
                    lcohs_trucking[i], generator_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], t_h2_storages[i] = plant.get_results(generators)
                    t_solar_capacities[i] = generator_capacities["Solar"]
                    t_wind_capacities[i] = generator_capacities["Wind"]
                else:
                    n = plant.get_hydrogen_network(hydrogen_demand_pipeline, wind_profile.time, country_series)
                    plant.set_wind_in_network(wind_potential, wind_max_capacity, country_series)
                    plant.set_pv_in_network(pv_potential, pv_max_capacity, country_series)  
                    plant.solve_model(solver)
                    lcohs_pipeline[i], generator_capacities, p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = plant.get_results(generators)
                    p_solar_capacities[i] = generator_capacities["Solar"]
                    p_wind_capacities[i] = generator_capacities["Wind"]
                
        # updating trucking hexagons
        hexagons[f'{demand_center} trucking solar capacity'] = t_solar_capacities
        hexagons[f'{demand_center} trucking wind capacity'] = t_wind_capacities
        hexagons[f'{demand_center} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{demand_center} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{demand_center} trucking H2 storage capacity'] = t_h2_storages
        # save trucking LCOH
        hexagons[f'{demand_center} trucking production cost'] = lcohs_trucking            

        # updating pipeline hexagons
        hexagons[f'{demand_center} pipeline solar capacity'] = p_solar_capacities
        hexagons[f'{demand_center} pipeline wind capacity'] = p_wind_capacities
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{demand_center} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{demand_center} pipeline H2 storage capacity'] = p_h2_storages
        # add optimal LCOH for each hexagon to hexagon file
        hexagons[f'{demand_center} pipeline production cost'] = lcohs_pipeline


    hexagons.to_file("results/hex.geojson", driver='GeoJSON', encoding='utf-8')