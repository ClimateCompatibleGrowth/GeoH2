from functions import CRF, cheapest_trucking_strategy, h2_conversion_stand, cheapest_pipeline_strategy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
from transport_optimization import calculate_road_construction_cost, calculate_distance_to_demand, check_folder_exists

if __name__ == "__main__":
    # ---------------------------------- Parameters variables ----------------------------------
    hexagons = gpd.read_file('data/hex_final_DJ.geojson') # SNAKEMAKE INPUT
    tech_params_filepath = 'parameters/technology_parameters.xlsx' # SNAKEMAKE INPUT
    demand_params_filepath = 'parameters/demand_parameters.xlsx' # SNAKEMAKE INPUT
    country_params_filepath = 'parameters/country_parameters.xlsx' # SNAKEMAKE INPUT
    conversion_params_filepath = 'parameters/conversion_parameters.xlsx' # SNAKEMAKE INPUT
    transport_params_filepath = 'parameters/transport_parameters.xlsx' # SNAKEMAKE INPUT
    pipeline_params_filepath = 'parameters/pipeline_parameters.xlsx' # SNAKEMAKE INPUT


    infra_data = pd.read_excel(tech_params_filepath,
                           sheet_name='Infra',
                           index_col='Infrastructure')
    
    demand_center_list = pd.read_excel(demand_params_filepath,
                                    index_col='Demand center',
                                    )
    demand_centers = demand_center_list.index
    water_data = pd.read_excel(tech_params_filepath,
                                sheet_name='Water',
                                index_col='Parameter'
                                ).squeeze("columns")
    country_params = pd.read_excel(country_params_filepath,
                                        index_col='Country')
    # ------------------------------------------------------------------------------------------

    # --------------------------------- Transport-optimization variables -----------------------
    needs_pipeline_construction = True # CONFIG YAML
    needs_road_construction = True # CONFIG YAML

    long_road_capex = infra_data.at['Long road','CAPEX']
    short_road_capex = infra_data.at['Short road','CAPEX']
    road_opex = infra_data.at['Short road','OPEX']
    # ------------------------------------------------------------------------------------------

    # ---------------------------------------- Water-cost variables ----------------------------------------
    h2o_costs_dom_water_bodies = np.empty(len(hexagons))
    h2o_costs_ocean = np.empty(len(hexagons))
    min_h2o_costs = np.empty(len(hexagons))

    electricity_demand_h2o_treatment = water_data['Freshwater treatment electricity demand (kWh/m3)']
    electricity_demand_ocean_h2o_treatment = water_data['Ocean water treatment electricity demand (kWh/m3)']
    water_transport_costs = water_data['Water transport cost (euros/100 km/m3)']
    water_spec_cost = water_data['Water specific cost (euros/m3)']
    water_demand = water_data['Water demand  (L/kg H2)']
    elec_price = country_params['Electricity price (euros/kWh)'].iloc[0]
    count = 0
    # ------------------------------------------------------------------------------------------------------

    check_folder_exists("results")
    # --------------------------------- Transport-optimization section ---------------------------------

    #%% calculate cost of hydrogen state conversion and transportation for demand
    # loop through all demand centers-- limit this on continential scale
    for demand_center in demand_centers:
        # Demand location based variables
        demand_center_lat = demand_center_list.loc[demand_center,'Lat [deg]']
        demand_center_lon = demand_center_list.loc[demand_center,'Lon [deg]']
        demand_location = Point(demand_center_lon, demand_center_lat)
        hydrogen_quantity = demand_center_list.loc[demand_center,'Annual demand [kg/a]']
        demand_state = demand_center_list.loc[demand_center,'Demand state']

        # Storage hexagons for costs calculated in the next for loop
        road_construction_costs = np.empty(len(hexagons))
        trucking_states = np.empty(len(hexagons),dtype='<U10')
        trucking_costs = np.empty(len(hexagons))
        pipeline_costs = np.empty(len(hexagons))

        # Prices from the country excel file
        elec_price = country_params['Electricity price (euros/kWh)'].iloc[0]
        heat_price = country_params['Heat price (euros/kWh)'].iloc[0]
        plant_interest_rate = country_params['Plant interest rate'].iloc[0]
        infrastructure_interest_rate = country_params['Infrastructure interest rate'].iloc[0]
        infrastructure_lifetime = country_params['Infrastructure lifetime (years)'].iloc[0]
        
        if demand_state not in ['500 bar','LH2','NH3']:
            raise NotImplementedError(f'{demand_state} demand not supported.')
        
        for i in range(len(hexagons)):
            distance_to_road = hexagons['road_dist'][i]
            hex_geometry = hexagons['geometry'][i]
            dist_to_demand = calculate_distance_to_demand(hex_geometry, demand_center_lat, demand_center_lon)

            #!!! maybe this is the place to set a restriction based on distance to demand center-- for all hexagons with a distance below some cutoff point
            # label demand location under consideration
            if hex_geometry.contains(demand_location) == True:
                # Calculate cost of converting hydrogen to a demand state for local demand (i.e. no transport)
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
                    continue
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
                    continue

            # determine elec_cost at demand to determine potential energy costs
            # Calculate cost of constructing a road to each hexagon
            if needs_road_construction == True:
                if distance_to_road==0:
                    road_construction_costs[i] = 0.
                elif distance_to_road!=0 and distance_to_road<10:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(distance_to_road, 
                                                         short_road_capex, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime, 
                                                         road_opex)
                else:
                    road_construction_costs[i] = \
                        calculate_road_construction_cost(distance_to_road, 
                                                         long_road_capex, 
                                                         infrastructure_interest_rate, 
                                                         infrastructure_lifetime, 
                                                         road_opex)
                    
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
            elif distance_to_road==0:
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
            elif distance_to_road>0: 
                trucking_costs[i]=trucking_states[i] = np.nan

            # Calculate costs of constructing a pipeline to each hexagon
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
            # --------------------------------------------------------------------------------------

            # --------------------------- Water-cost section ---------------------------
            # Water cost for each hexagon for each kg hydrogen produced
            if count == 0:
                waterbody_dist = hexagons['waterbody_dist'][i]
                waterway_dist = hexagons['waterway_dist'][i]
                ocean_dist = hexagons['ocean_dist'][i]
                
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
            
            # --------------------------------------------------------------------------
        count+=1
        # ---------------------------- Updating transport-optimization section ----------------------------
        
        # Hexagon file updated with each demand center's costs and states
        hexagons[f'{demand_center} road construction costs'] = road_construction_costs/hydrogen_quantity
        hexagons[f'{demand_center} trucking transport and conversion costs'] = trucking_costs # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} trucking state'] = trucking_states # cost of road construction, supply conversion, trucking transport, and demand conversion
        hexagons[f'{demand_center} pipeline transport and conversion costs'] = pipeline_costs # cost of supply conversion, pipeline transport, and demand conversion
        # --------------------------------------------------------------------------------------------------
    # ---------------------- Updating water-cost section ---------------------
    hexagons['Ocean water costs'] = h2o_costs_ocean
    hexagons['Freshwater costs'] = h2o_costs_dom_water_bodies
    hexagons['Lowest water cost'] = min_h2o_costs
    # ------------------------------------------------------------------------

    # ------------------------ Updating total-costs section ------------------------
    for demand_center in demand_centers:
        hexagons[f'{demand_center} trucking total cost'] =\
            hexagons[f'{demand_center} road construction costs'] +\
                hexagons[f'{demand_center} trucking transport and conversion costs'] +\
                    hexagons['Lowest water cost']
                        # hexagons[f'{demand_center} trucking production cost'] +\ HYRDOGEN OPTIMIZATION
        hexagons[f'{demand_center} pipeline total cost'] =\
                hexagons[f'{demand_center} pipeline transport and conversion costs'] +\
                    hexagons['Lowest water cost']
                        # hexagons[f'{demand_center} pipeline production cost'] +\ HYDROGEN OPTIMIZATION

                        
        for i in range(len(hexagons)):
            hexagons.loc[i, f'{demand_center} lowest cost'] = np.nanmin(
                                    [hexagons.loc[i, f'{demand_center} trucking total cost'],
                                    hexagons.loc[i, f'{demand_center} pipeline total cost']
                                    ])
    # ------------------------------------------------------------------------------

    hexagons.to_file("results/completed_hex_DJ.geojson", driver='GeoJSON', encoding='utf-8')