import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from functions import cheapest_trucking_strategy, h2_conversion_stand, cheapest_pipeline_strategy
from network import Network
from plant_optimization import get_demand_schedule, get_generator_profile, solve_model, get_results
from transport_optimization import calculate_road_construction_cost, calculate_dist_to_demand
from utils import check_folder_exists

if __name__ == "__main__":
    # ---------------------------------- Parameters variables ----------------------------------
    hexagons = gpd.read_file('data/hex_final_DJ.geojson') # SNAKEMAKE INPUT
    tech_params_filepath = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    conversion_params_filepath = 'parameters/conversion_parameters.xlsx' # SNAKEMAKE INPUT
    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    pipeline_params_filepath = 'parameters/pipeline_parameters.xlsx' # SNAKEMAKE INPUT

    # This is probably just me but I find it easier to read on one line than spaced down (unless too long)
    infra_data = pd.read_excel(tech_params_filepath, sheet_name='Infra', index_col='Infrastructure')
    demand_params = pd.read_excel(demand_params_filepath,index_col='Demand center',)
    demand_centers = demand_params.index
    water_data = pd.read_excel(tech_params_filepath, sheet_name='Water', index_col='Parameter').squeeze("columns")
    country_params = pd.read_excel(country_params_filepath, index_col='Country')
    elec_price = country_params['Electricity price (euros/kWh)'].iloc[0]

    len_hexagons = len(hexagons)

    # --------------------------------- Transport-optimization variables -----------------------
    needs_pipeline_construction = True # CONFIG YAML
    needs_road_construction = True # CONFIG YAML

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']

    # ---------------------------------------- Water-cost variables ----------------------------------------
    h2o_costs_dom_water_bodies = np.empty(len_hexagons)
    h2o_costs_ocean = np.empty(len_hexagons)
    min_h2o_costs = np.empty(len_hexagons)

    electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
    electricity_demand_ocean_h2o_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
    water_transport_costs = water_data['Water transport cost (euros/100 km/m3)']
    water_spec_cost = water_data['Water specific cost (euros/m3)']
    water_demand = water_data['Water demand  (L/kg H2)']
    demand_case_count = 0

    # ------------------------------------ Plant-optimization variables ------------------------------------
    country_series = country_params.iloc[0]
    weather_year = 2022             # snakemake.wildcards.weather_year
    end_weather_year = 2023         # int(snakemake.wildcards.weather_year)+1
    start_date = f'{weather_year}-01-01'
    end_date = f'{end_weather_year}-01-01'
    solver = "gurobi" # maybe make this into a snakemake wildcard? I think that makes the most sense!!
    generators = {"Wind" : [], "Solar" : []} # already in the config as list, used in map_costs.py line 268
    pipeline_construction = True # snakemake config

    cutout = atlite.Cutout('cutouts/DJ_2022.nc') # SNAKEMAKE INPUT
    # Get a uniform capacity layout for all grid cells. https://atlite.readthedocs.io/en/master/ref_api.html
    # Alycia to double-check we are using the right layout
    layout = cutout.uniform_layout()
    profiles = []

    for gen in generators.keys():
         profiles.append(get_generator_profile(gen, cutout, layout, hexagons))

    times = profiles[0].time
    # ---------------------------------------- Execution starts here --------------------------------------------

    check_folder_exists("results")

    # Calculate cost of hydrogen state conversion and transportation for demand

    # Loop through all demand centers-- limit this on continental scale
    for demand_center in demand_centers:
        # ------------------------------ Transport-optimization section part 1 ------------------------------
        # Get demand-location-based variables
        demand_center_lat = demand_params.loc[demand_center,'Lat [deg]']
        demand_center_lon = demand_params.loc[demand_center,'Lon [deg]']
        demand_location = Point(demand_center_lon, demand_center_lat)
        hydrogen_quantity = demand_params.loc[demand_center,'Annual demand [kg/a]']
        demand_state = demand_params.loc[demand_center,'Demand state']

        # Storage hexagons for costs calculated in the next for loop
        road_construction_costs = np.empty(len_hexagons)
        trucking_states = np.empty(len_hexagons,dtype='<U10')
        trucking_costs = np.empty(len_hexagons)
        pipeline_costs = np.empty(len_hexagons)

        # Does this need to be inside the demand centre loop? Or could they move out? YES - has been moved out in its separate file
        # Get prices from the country excel file
        heat_price = country_params['Heat price (euros/kWh)'].iloc[0]
        plant_interest_rate = country_params['Plant interest rate'].iloc[0]
        infrastructure_interest_rate = country_params['Infrastructure interest rate'].iloc[0]
        infrastructure_lifetime = country_params['Infrastructure lifetime (years)'].iloc[0]
        
        if demand_state not in ['500 bar','LH2','NH3']:
            raise NotImplementedError(f'{demand_state} demand not supported.')

        # --------------------------------- Plant-optimization section  part 1 ---------------------------------
        # Make places to store trucking results
        lcohs_trucking = np.zeros(len_hexagons)
        t_generators_capacities = { "Wind" : [], "Solar" : []} # snakemake config
        t_electrolyzer_capacities= np.zeros(len_hexagons)
        t_battery_capacities = np.zeros(len_hexagons)
        t_h2_storages= np.zeros(len_hexagons)
        
        # Make places to store pipeline results
        lcohs_pipeline = np.zeros(len_hexagons)
        p_generators_capacities = { "Wind" : [], "Solar" : []} # snakemake config
        p_electrolyzer_capacities= np.zeros(len_hexagons)
        p_battery_capacities = np.zeros(len_hexagons)
        p_h2_storages= np.zeros(len_hexagons)

        # Get the quantity of demand at the demand center
        hydrogen_quantity = demand_params.loc[demand_center,'Annual demand [kg/a]']

        # Loop through all hexagons
        for i in range(len_hexagons):
            # ------------------------------ Transport-optimization section part 2 ------------------------------
            # Get the hexagon-specific parameters
            distance_to_road = hexagons['road_dist'][i]
            hex_geometry = hexagons['geometry'][i]
            dist_to_demand = calculate_dist_to_demand(hex_geometry, demand_center_lat, demand_center_lon)

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration

            # If the hexagon contains the demand location
            if hex_geometry.contains(demand_location) == True:
                # The trucking cost is just the conversion cost for local use
                if demand_state == 'NH3':
                    trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state+'_load',
                                             hydrogen_quantity,
                                             elec_price,
                                             heat_price,
                                             plant_interest_rate,
                                             conversion_params_filepath
                                             )[2]/hydrogen_quantity
                    trucking_states[i] = "None"
                    road_construction_costs[i] = 0.
                else:
                    trucking_costs[i]=pipeline_costs[i]=h2_conversion_stand(demand_state,
                                             hydrogen_quantity,
                                             elec_price,
                                             heat_price,
                                             plant_interest_rate,
                                             conversion_params_filepath
                                             )[2]/hydrogen_quantity
                    trucking_states[i] = "None"
                    road_construction_costs[i] = 0.

            # Otherwise, if the hexagon does not contain the demand center
            # This structure is already in the transport file - move it identically over here.
            else:
            # Calculate the cost of constructing a road to the hexagon if needed
                if needs_road_construction == True:
                    # If there 0 distance to road, there is no road construction cost.
                    if distance_to_road==0:
                        road_construction_costs[i] = 0.
                    # If the distance is more than 0 and less than 10, use short road costs.
                    elif distance_to_road!=0 and distance_to_road<10:
                        road_construction_costs[i] = \
                            calculate_road_construction_cost(distance_to_road,
                                                             short_road_capex,
                                                             infrastructure_interest_rate,
                                                             infrastructure_lifetime,
                                                             road_opex)
                    else:   # Otherwise (i.e., if distance is more than 10), use long road costs.
                        road_construction_costs[i] = \
                            calculate_road_construction_cost(distance_to_road,
                                                             long_road_capex,
                                                             infrastructure_interest_rate,
                                                             infrastructure_lifetime,
                                                             road_opex)
                    # Then find cheapest trucking strategy.
                    trucking_costs[i], trucking_states[i] =\
                        cheapest_trucking_strategy(demand_state,
                                                    hydrogen_quantity,
                                                    dist_to_demand,
                                                    elec_price,
                                                    heat_price,
                                                    infrastructure_interest_rate,
                                                    conversion_params_filepath,
                                                    transport_params_filepath,
                                                    )
                # Otherwise, if road construction not allowed:
                else:
                    # If distance to road is 0, just get cheapest trucking strategy
                    if distance_to_road==0:
                        trucking_costs[i], trucking_states[i] =\
                            cheapest_trucking_strategy(demand_state,
                                                        hydrogen_quantity,
                                                        dist_to_demand,
                                                        elec_price,
                                                        heat_price,
                                                        infrastructure_interest_rate,
                                                        conversion_params_filepath,
                                                        transport_params_filepath,
                                                        )
                    # And if road construction is not allowed and distance to road is > 0, trucking states are nan.
                    # Sam to confirm whether assigning nan will cause future issues in code.
                    elif distance_to_road>0:
                        trucking_costs[i]=trucking_states[i] = np.nan

                # Calculate costs of constructing a pipeline to the hexagon if allowed
                if needs_pipeline_construction== True:
                    pipeline_cost, pipeline_type =\
                        cheapest_pipeline_strategy(demand_state,
                                                hydrogen_quantity,
                                                dist_to_demand,
                                                elec_price,
                                                heat_price,
                                                infrastructure_interest_rate,
                                                conversion_params_filepath,
                                                pipeline_params_filepath,
                                                )
                    pipeline_costs[i] = pipeline_cost
                else:
                    pipeline_costs[i] = np.nan

            # --------------------------------- Plant-optimization section  part 2 ---------------------------------
            trucking_state = trucking_states[i]
            gen_index = 0
            # Is there no zero-indexed hexagon? Ahh - or just reset after the first one?
            if i > 0:
                generators = { "Wind" : [], "Solar" : []} # snakemake config
            # Get the max capacity for each generation type
            for gen in generators.keys():
                potential = profiles[gen_index].sel(hexagon = i)
                # Eventually make a for loop - we can change the theo_turbines name to be Wind
                if gen == "Wind":
                    max_capacity = hexagons.loc[i,'theo_turbines']*4 # We'll need to remove this hard-coded 4 eventually CONFIG FILE - 4 MW turbine in spatial data prep
                elif gen == "Solar":
                    max_capacity = hexagons.loc[i,'theo_pv']
                # Eventually move loops to something like this so we don't have ifs - max_capacity = hexagons.loc[i, gen] * SNAKEMAKE_CONFIG_GEN_SIZE

                generators[gen].append(potential)
                generators[gen].append(max_capacity)
                gen_index += 1

            # Get the demand schedule for both pipeline and trucking transport
            trucking_demand_schedule, pipeline_demand_schedule =\
                get_demand_schedule(hydrogen_quantity,
                                start_date,
                                end_date,
                                trucking_state,
                                transport_params_filepath)

            # For each transport type, set up the network
            transport_types = ["trucking", "pipeline"]
            for transport in transport_types:
                network = Network("Hydrogen", generators)
                # For trucking, set it up with the trucking demand schedule
                if transport == "trucking":
                    if trucking_state != np.nan: # For all where there is a trucking state that's viable
                        network.set_network(trucking_demand_schedule, times, country_series)
                        network.set_generators_in_network(country_series)
                        solve_model(network.n, solver)
                        lcohs_trucking[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], t_h2_storages[i] = get_results(network.n, trucking_demand_schedule, generators)
                        for gen, capacity in generators_capacities.items():
                            t_generators_capacities[gen].append(capacity)
                    else: # If the hexagon has no viable trucking state (i.e., no roads reach it), set everything to nan.
                        lcohs_trucking[i], generators_capacities, t_electrolyzer_capacities[i], t_battery_capacities[i], \
                        t_h2_storages[i] = np.nan
                        for gen, capacity in generators_capacities.items():
                            t_generators_capacities[gen].append(np.nan)
                # For pipeline, set it up with pipeline demand schedule if construction is true
                else:
                    if pipeline_construction == True:
                        network.set_network(pipeline_demand_schedule, times, country_series)
                        network.set_generators_in_network(country_series)
                        solve_model(network.n, solver)
                        lcohs_pipeline[i], generators_capacities, p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = get_results(network.n, pipeline_demand_schedule, generators)
                        for gen, capacity in generators_capacities.items():
                            p_generators_capacities[gen].append(capacity)
                    # If construction is false, you can't transport it - everything gets nan UNLESS you're in the demand centre hexagon
                    else:
                        if hex_geometry.contains(demand_location) == True:
                            network.set_network(pipeline_demand_schedule, times, country_series)
                            network.set_generators_in_network(country_series)
                            solve_model(network.n, solver)
                            lcohs_pipeline[i], generators_capacities, p_electrolyzer_capacities[i], \
                            p_battery_capacities[i], p_h2_storages[i] = get_results(network.n, pipeline_demand_schedule,
                                                                                    generators)
                            for gen, capacity in generators_capacities.items():
                                p_generators_capacities[gen].append(capacity)
                        else:
                            lcohs_pipeline[i], p_electrolyzer_capacities[i], p_battery_capacities[i], p_h2_storages[i] = np.nan
                            for gen in p_generators_capacities.keys():
                                p_generators_capacities[gen].append(np.nan)

            # --------------------------- Water-cost section ---------------------------
            # Water cost for each hexagon for each kg hydrogen produced
            # Water cost is worked out once for each hexagon
            if demand_case_count == 0:
                waterbody_dist = hexagons['waterbody_dist'][i]
                waterway_dist = hexagons['waterway_dist'][i]
                ocean_dist = hexagons['ocean_dist'][i]
                # I assume the 100 and 1000 are just units?
                h2o_costs_dom_water_bodies[i] =(water_spec_cost 
                                                    + (water_transport_costs/100)
                                                    * min(waterbody_dist, waterway_dist) 
                                                    + electricity_demand_h2o_treatment
                                                    * elec_price
                                                    ) * water_demand/1000
                h2o_costs_ocean[i] =(water_spec_cost 
                                        + (water_transport_costs/100)
                                        * ocean_dist
                                        + electricity_demand_ocean_h2o_treatment
                                        * elec_price
                                        ) * water_demand/1000
                
                min_h2o_costs[i] = min(h2o_costs_dom_water_bodies[i], h2o_costs_ocean[i])

        demand_case_count+=1
        # ---------------------------- Updating hexagon with transport results ----------------------------
        
        # Updating hexagon file updated with each demand center's costs and states
        hexagons[f'{demand_center} road construction costs'] = road_construction_costs/hydrogen_quantity
        hexagons[f'{demand_center} trucking transport and conversion costs'] = trucking_costs # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} trucking state'] = trucking_states # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} pipeline transport and conversion costs'] = pipeline_costs # cost of supply conversion, pipeline transport, and demand conversion

        # Updating trucking-based results in hexagon file
        for gen, capacities in t_generators_capacities.items():
            hexagons[f'{demand_center} trucking {gen.lower()} capacity'] = capacities
        hexagons[f'{demand_center} trucking electrolyzer capacity'] = t_electrolyzer_capacities
        hexagons[f'{demand_center} trucking battery capacity'] = t_battery_capacities
        hexagons[f'{demand_center} trucking H2 storage capacity'] = t_h2_storages
        hexagons[f'{demand_center} trucking production cost'] = lcohs_trucking            

        # Updating pipeline-based results in hexagon file
        for gen, capacities in p_generators_capacities.items():
            hexagons[f'{demand_center} pipeline {gen.lower()} capacity'] = capacities
        hexagons[f'{demand_center} pipeline electrolyzer capacity'] = p_electrolyzer_capacities
        hexagons[f'{demand_center} pipeline battery capacity'] = p_battery_capacities
        hexagons[f'{demand_center} pipeline H2 storage capacity'] = p_h2_storages
        hexagons[f'{demand_center} pipeline production cost'] = lcohs_pipeline

    # ---------------------- Updating water results in hexagon file ---------------------
    hexagons['Ocean water costs'] = h2o_costs_ocean
    hexagons['Freshwater costs'] = h2o_costs_dom_water_bodies
    hexagons['Lowest water cost'] = min_h2o_costs

    # ------------------------ Updating total cost results in hexagon file ------------------------
    # Get lowest cost for each transport type
    for demand_center in demand_centers:
        hexagons[f'{demand_center} trucking total cost'] =\
            hexagons[f'{demand_center} road construction costs'] +\
                hexagons[f'{demand_center} trucking transport and conversion costs'] +\
                        hexagons[f'{demand_center} trucking production cost'] +\
                            hexagons['Lowest water cost']
        hexagons[f'{demand_center} pipeline total cost'] =\
                hexagons[f'{demand_center} pipeline transport and conversion costs'] +\
                    hexagons[f'{demand_center} pipeline production cost'] +\
                        hexagons['Lowest water cost']
        # Get the lowest between the trucking and the pipeline options
        for i in range(len_hexagons):
            hexagons.loc[i, f'{demand_center} lowest cost'] = np.nanmin(
                                    [hexagons.loc[i, f'{demand_center} trucking total cost'],
                                    hexagons.loc[i, f'{demand_center} pipeline total cost']
                                    ])
    # Save all results to file
    hexagons.to_file("results/completed_hex_DJ.geojson", driver='GeoJSON', encoding='utf-8')